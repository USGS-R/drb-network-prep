get_nhdv2_flowlines <- function(huc8){
  #' 
  #' @description Function to download NHDPlusV2 flowlines for specified HUC8 watersheds
  #'
  #' @param huc8 string or character vector containing the HUC8 watersheds of interest
  #'
  #' @value An sf object containing NHDPlusV2 flowlines covering the HUC8 watersheds of interest
  
  # Download flowlines for each HUC8
  flines_by_huc <- huc8 %>%
    lapply(.,function(x){
      # Create spatial object of huc8 basin
      huc8_basin <- suppressMessages(nhdplusTools::get_huc8(id=x))
      # Create spatial object of huc8 bbox
      huc8_bbox <- sf::st_as_sfc(st_bbox(huc8_basin))
      # Download NHDPlusV2 flowlines within each huc8 bbox
      huc8_flines <- suppressMessages(nhdplusTools::get_nhdplus(AOI=huc8_bbox,
                                                                realization="flowline"))
    })
  
  # Bind HUC8 flowlines together 
  flines <- flines_by_huc %>%
    bind_rows() %>%
    # Reformat variable names to uppercase
    rename_with(.,toupper,id:enabled) %>%
    # Remove duplicated COMIDs that result from retrieving flowlines associated 
    # with (overlapping) HUC8 bbox's
    group_by(COMID) %>%
    slice(1) %>%
    ungroup() %>%
    # create a new column that uses the REACHCODE attribute to assign
    # which huc8 sub-basin each flowline reach belongs to
    mutate(huc8_code = str_sub(REACHCODE, 0, 8))
  
  # Filter flowlines to only include reaches within the basin/watershed of 
  # interest (as defined by the huc8 codes provided)
  flines_hucs <- flines %>%
    filter(huc8_code %in% huc8)
  
  # Check whether any other flowlines drain to the target huc8 sub-basins,
  # even if the reachcode suggests otherwise
  flines_others <- flines %>%
    filter(!huc8_code %in% huc8,
           LEVELPATHI %in% flines_hucs$LEVELPATHI)
  
  # Return flowlines
  flines_out <- bind_rows(flines_hucs, flines_others)
    
  return(flines_out)
  
}



#' Function to download NHDPlusv2 catchments given a set of COMIDs
#' 
#' @param comid character string or vector of character strings containing
#' the common identifier (COMID or featureid) of the desired catchment(s).
#' 
get_nhdplusv2_catchments <- function(comid){
  
  # Chunk desired COMIDs into groups, where each group has no more than
  # 50 COMID's to avoid timeout errors when downloading nhdplus subsets
  # using helper functions from nhdplusTools
  comid_df <- tibble(COMID = comid) %>%
    mutate(comid_n = row_number(),
           download_grp = ((comid_n -1) %/% 50) + 1)
  
  # Download catchments associated with desired COMIDs and return an
  # sf data frame containing the polygons and value-added attributes.
  # suppress messages "Found invalid geometry, attempting to fix" from
  # printing to the console.
  catchments <- comid_df %>%
    split(., .$download_grp) %>%
    lapply(., function(x){
      cats_sub <- suppressMessages(nhdplusTools::get_nhdplus(comid = x$COMID, 
                                                             realization = "catchment"))
    }) %>%
    bind_rows()
  
  # Format variable names
  catchments_out <- catchments %>%
    rename_with(., toupper, id:shape_area) %>%
    rename(COMID = FEATUREID)
  
  # Inform user how many catchments were returned
  message(sprintf("Returning %s catchments of %s total COMID's requested",
                  format(length(catchments_out$ID), big.mark=","), 
                  format(length(comid), big.mark=",")))
  
  return(catchments_out)
  
}

