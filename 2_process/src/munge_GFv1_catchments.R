find_intersecting_hru <- function(prms_line,prms_hrus){
  #' 
  #' @description This function finds all of the HRU polygons that intersect a
  #' PRMS segment, and returns the HRU ID (i.e., hru_segment) that has the
  #' greatest overlap with the PRMS segment
  #' 
  #' @param prms_line sf object containing one PRMS river segment of interest
  #' prms_line must contain variables subsegid and geometry
  #' @param prms_hrus sf (multi)polygon containing the HRU's from the NHGFv1
  #
  
  # Transform PRMS segment and HRU's into a consistent projection
  # using Albers Equal Area Conic CONUS
  prms_line_proj <- sf::st_transform(prms_line, 5070) %>%
    sf::st_buffer(5) %>%
    select(subsegid, subsegseg, subseglen, fromsegs, toseg, tosubseg, segidnat)
  prms_hrus_proj <- sf::st_transform(prms_hrus, 5070)

  # Return the intersection between prms_line and HRU polygons 
  # suppress warnings that "attribute variables are assumed to be spatially constant
  # throughout all geometries"
  prms_line_intsct_hrus <- sf::st_intersection(prms_line_proj,prms_hrus_proj) %>%
    suppressWarnings()
  
  # Find which hru_segment maximizes the area of overlap between the HRU and prms_line
  prms_line_intsct_hrus$area_overlap <- as.numeric(sf::st_area(prms_line_intsct_hrus))
  
  # Return attribute `hru_segment` for the HRU that has the greatest overlap with prms_line
  hru_segment <- prms_line_intsct_hrus$hru_segment[which.max(prms_line_intsct_hrus$area_overlap)]
  
  return(hru_segment)
  
}



return_nhdv2_cat_shps <- function(prms_line,segs_w_comids,verbose = TRUE){
  #' 
  #' @description This function takes a data table of PRMS segment ID's and
  #' its contributing NHDv2 tributary catchments and returns those catchments
  #' as an sf object
  #' 
  #' @param prms_line sf object containing one PRMS river segment of interest
  #' prms_line must contain variables subsegid and geometry
  #' @param segs_w_comids data frame containing the PRMS segment ids and the contributing COMID's
  #' segs_w_comids must contain variables PRMS_segid and COMID
  #' @param verbose logical, defaults to TRUE; if FALSE do not include diagnostic messages
  #' 
  #' @details sf > v1.0 causes issues with some NHD polygons and spherical coordinates, see
  #' info on the change: https://r-spatial.org/r/2020/06/17/s2.html#sf-10-goodbye-flat-earth-welcome-s2-spherical-geometry
  #' 
  #' To avoid issues with NHD, we won't use the S2 engine for prepping NHD catchments. This
  #' function toggles off the use of S2 and then resets the global options when the function is
  #' complete.
  #' 
  
  sf::sf_use_s2(FALSE)
  on.exit(sf::sf_use_s2(TRUE))
  
  # Subset segs_w_comids to only include the PRMS segment of interest
  nhd_cats <- segs_w_comids %>%
    filter(PRMS_segid == prms_line$subsegid)
  
  # Use nhdplusTools to retrieve the catchment polygons for COMIDs of interest
  nhd_cats_sf <- nhdplusTools::get_nhdplus(comid = c(nhd_cats$COMID),realization = "catchment",t_srs = 5070) 
  
  if(verbose == TRUE){
    message(sprintf("Gathering NHDv2 catchments to create HRUs for subsegid %s; %s COMIDs expected and %s COMIDs returned",
                    unique(prms_line$subsegid),length(nhd_cats$COMID),length(nhd_cats_sf$featureid)))
  }
  
  # Format the catchment polygons
  nhd_cats_sf_out <- nhd_cats_sf %>%
    # dissolve boundaries between individual NHDv2 catchments
    sf::st_union() %>%
    # convert back to sf data.frame object and add the segment ID
    sf::st_as_sf() %>%
    mutate(PRMS_segid = unique(prms_line$subsegid)) %>%
    rename(geometry = x)

  return(nhd_cats_sf_out)
  
}



munge_GFv1_catchments <- function(prms_lines,prms_hrus,segs_w_comids,segs_split,crs_out = 4326,verbose = TRUE){
  #' 
  #' @description This function munges the PRMS HRU's (~ the PRMS catchments) so that every
  #' PRMS segment in the Delaware River Basin has at least one corresponding HRU and so that 
  #' the catchment boundaries are adjusted for NHGFv1 segments that were split in the 
  #' `delaware-model-prep` pipeline (https://github.com/USGS-R/delaware-model-prep) 
  #' 
  #' @param prms_lines sf object containing the PRMS river segments for the DRB
  #' prms_lines must contain variables subsegid and geometry
  #' @param segs_w_comids data frame containing the PRMS segment ids and the contributing COMID's
  #' segs_w_comids must contain variables PRMS_segid and COMID
  #' @param prms_hrus sf (multi)polygon containing the HRU's from the NHGFv1
  #' @param segs_split character vector indicating the PRMS_segid of the segments that were split in 
  #' the delaware-model-prep pipeline
  #' @param crs_out numeric EPSG code indicating desired crs. Defaults to EPSG: 4326, WGS84.
  #' @param verbose logical, defaults to TRUE; if FALSE do not include diagnostic messages
  #
  
  # 1. Fill PRMS lines data frame so that every segment has a value for the attribute
  # `hru_segment`, which can be used as a key to join each segment to its corresponding HRU polygons.
  prms_hrus_filled <- prms_lines %>%
    # Transform prms_lines into a list with a data frame for each PRMS segment (represented by a unique subsegid) 
    split(.,.$subsegid) %>%
    # For each PRMS segment, find the HRU with the greatest overlap and add as an attribute
    lapply(.,function(x){
      x$hru_segment <- ifelse(x$subsegseg %in% prms_hrus$hru_segment,
                              x$subsegseg,
                              find_intersecting_hru(x,prms_hrus))
      return(x)
    }) %>%
    # Reformat list into a data frame with one row per PRMS segment
    bind_rows()
  
  # 2. For select segments, redefine the catchment polygon using the contributing NHDv2 catchments
  prms_splitsegs <- prms_hrus_filled %>%
    # Filter PRMS segments to only include select segments that require special handling
    filter(subsegid %in% segs_split) %>%
    # Transform PRMS segments into a list with a data frame for each segment
    split(.,.$subsegid) %>%
    # For each segment requiring special handling, return the polygon that represents the catchment area
    lapply(.,return_nhdv2_cat_shps, segs_w_comids = segs_w_comids, verbose = verbose) %>%
    # Reformat list into a data frame with one row per PRMS segment
    bind_rows()
  
  # 3. Concatenate munged catchment polygons
  catchments_ls <- list()
  
  # Fill empty list with appropriate catchments
  for(i in seq_along(prms_hrus_filled$subsegid)){

    prms_line <- prms_hrus_filled[i,]
    
    # If PRMS segment is one of the segments that require special handling, 
    # retrieve NHDPlusV2 catchments and use those as the catchment boundary.
    # Else, represent the catchment boundaries using intersecting HRU polygon.
    if(prms_line$subsegid %in% segs_split){
      cat <- prms_splitsegs %>%
        filter(PRMS_segid == prms_line$subsegid) %>%
        # used manual inspection to match updated catchments with closest
        # hru_id_nat from the NHM v1.0
        mutate(hru_segment = prms_line$subsegid,
               hru_id_reg = case_when(
                 PRMS_segid == "3_1" ~ "4042",
                 PRMS_segid == "3_2" ~ "4043",
                 PRMS_segid == "8_1" ~ "4785",
                 PRMS_segid == "8_2" ~ "4785",
                 PRMS_segid == "51_1" ~ "3433",
                 PRMS_segid == "51_2" ~ "3428"),
               hru_id_nat = case_when(
                 PRMS_segid == "3_1" ~ "6504",
                 PRMS_segid == "3_2" ~ "6505",
                 PRMS_segid == "8_1" ~ "7247",
                 PRMS_segid == "8_2" ~ "7247",
                 PRMS_segid == "51_1" ~ "5895",
                 PRMS_segid == "51_2" ~ "5890"),
               region = "02",
               hru_area_m2 = as.numeric(sf::st_area(.)),
               hru_area_km2 = round(hru_area_m2/(1000*1000),3)) %>%
        select(PRMS_segid,hru_segment,hru_id_reg,hru_id_nat,hru_area_m2,hru_area_km2,region)
    } else {
      cat <- prms_hrus %>%
        filter(hru_segment == prms_line$hru_segment) %>%
        mutate(PRMS_segid = prms_line$subsegid,
               hru_area_km2 = round(Shape_Area/(1000*1000),3),
               hru_id_reg = as.character(hru_id_reg),
               hru_id_nat = as.character(hru_id_nat),
               hru_segment = as.character(hru_segment)) %>%
        select(c(PRMS_segid,hru_id_reg,Shape_Area,hru_area_km2,
                 any_of(gf_cols_select))) %>%
        relocate(c(hru_segment,hru_id_nat),.after = PRMS_segid) %>%
        rename(geometry = Shape, hru_area_m2 = Shape_Area) 
    }
    
    # Transform edited catchments to user-specified coordinate reference system
    cat_out <- cat %>%
      sf::st_transform(.,crs_out) 
    
    # Append edited catchments
    catchments_ls[[i]] <- cat_out
    
  }
  
  # Bind edited catchments into a single data frame
  catchments_proc <- bind_rows(catchments_ls)

  return(catchments_proc)
  
}

