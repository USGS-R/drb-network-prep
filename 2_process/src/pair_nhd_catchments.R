pair_nhd_catchments <- function(prms_line, prms_hrus, min_area_overlap, xwalk_table,
                                nhd_lines, drb_segs_spatial, omit_zero_area_flines,
                                zero_area_flines_max_length){
  #'
  #' @description This function finds the NHDPlusV2 catchments that overlap the HRU's for each PRMS segment in the DRB. 
  #' For the majority of DRB segments, this function finds the catchments that directly drain to a given segment 
  #' by navigating the NHDPlusV2 value-added attributes tables.  
  #' For a subset of segments, the contributing area by NHD navigation differed considerably from the 
  #' contributing area represented by the corresponding HRU polygons (based on visual inspection). This
  #' function applies special handling for those select segments and uses a spatial join to 
  #' find the NHDPlusV2 catchments that overlap the PRMS HRU's.
  #' 
  #' @param prms_line sf linestring containing a target PRMS segment
  #' @param prms_hrus sf (multi)polygon containing the HRU's from the NHGF
  #' @param min_area_overlap float; value indicating the minimum proportion of NHDPlusV2 catchment area that 
  #' overlaps the PRMS polygon in order to be retained
  #' @param xwalk_table data frame containing the NHDPlusV2/PRMS reach-to-segment crosswalk information
  #' @param nhd_lines sf object containing NHDPlusV2 flowlines for area of interest
  #' nhd_lines must contain variables COMID, PATHLENGTH, LENGTHKM, and HYDROSEQ
  #' @param drb_segs_spatial character vector containing the identity of PRMS segments that require special handling
  #' For these segments, the contributing NHD catchments will be determined using a spatial join with the 
  #' PRMS HRU polygons
  #' @param omit_zero_area_flines logical; if TRUE, only return NHDPlusv2 reaches where the NHDPlus atttribute
  #' "AREASQKM" is greater than zero.
  #' @param zero_area_flines_max_length integer; zero-area reaches are expected to be
  #' relatively short, so a message is printed to the console if `omit_zero_area_flines` 
  #' is TRUE and `LENGTHKM` of an omitted reach exceeds `zero_area_flines_max_length`.
  #' 
  #' @examples pair_nhd_catchments(prms_line = p1_reaches_sf[p1_reaches_sf$subsegid == "3_1",],
  #'                               prms_hrus = p1_catchments_sf,
  #'                               min_area_overlap = 0.7,
  #'                               xwalk_table = p2_prms_nhdv2_xwalk,
  #'                               nhd_lines = p1_nhdv2reaches_sf,
  #'                               drb_segs_spatial = "3_1",
  #'                               omit_zero_area_flines = TRUE, 
  #'                               zero_area_flines_max_length = 0.25)
  #'                               
  if(prms_line$subsegid %in% drb_segs_spatial){
    # If PRMS segment requires special handling, find contributing NHDPlusV2 catchments
    # through a spatial join with the corresponding PRMS HRU's
    comid_cat <- get_intersecting_nhdplus_catchments(prms_line,
                                                     prms_hrus,
                                                     min_area_overlap, 
                                                     omit_zero_area_flines,
                                                     zero_area_flines_max_length)    
  } else {
    # Otherwise, find contributing NHDPlusV2 catchments by navigating upstream tributaries
    comid_cat <- get_upstream_nhd_catchments(nhd_lines,
                                             xwalk_table,
                                             prms_line,
                                             omit_zero_area_flines,
                                             zero_area_flines_max_length)
  }
  
  return(comid_cat)
  
}



get_intersecting_nhdplus_catchments <- function(prms_line, prms_hrus,
                                                min_area_overlap,
                                                omit_zero_area_flines,
                                                zero_area_flines_max_length){
  #' @description This function uses a spatial join to find the NHDPlusV2 catchments 
  #' that overlap the PRMS HRU's.
  #' 
  #' @param prms_line sf linestring containing a target PRMS segment
  #' @param prms_hrus sf (multi)polygon containing the HRU's from the NHGF
  #' @param min_area_overlap float; value indicating the minimum proportion of NHDPlusV2 
  #' catchment area that overlaps the PRMS polygon in order to be retained
  #' @param omit_zero_area_flines logical; if TRUE, only return NHDPlusv2 reaches where 
  #' the NHDPlus atttribute "AREASQKM" is greater than zero.
  #' @param zero_area_flines_max_length integer; zero-area reaches are expected to be
  #' relatively short, so a message is printed to the console if `omit_zero_area_flines` 
  #' is TRUE and `LENGTHKM` of an omitted reach exceeds `zero_area_flines_max_length`.
  #' 
  
  # Select corresponding HRU polygon 
  prms_poly <- prms_hrus %>% 
    filter(hru_segment == prms_line$subsegseg) %>% 
    st_union()
  
  if(length(st_geometry(prms_poly))>0){
    
    # Download NHDPlusV2 catchments and calculate original area
    nhd_cats <- nhdplusTools::get_nhdplus(AOI = prms_poly,realization = "catchment") %>%
      suppressMessages() %>%
      mutate(area_orig = as.numeric(st_area(.)))
    
    # Find intersection between NHD catchments and HRU polygon and calculate *intersecting* area
    nhd_cats_intrsct <- st_intersection(nhd_cats,prms_poly) %>%
      suppressWarnings() %>%
      mutate(area_int = as.numeric(st_area(.)),
             area_prop = as.numeric(area_int/area_orig)) %>%
      # retain NHD catchments that overlap >= X% of their area with the HRU polygon
      filter(area_prop >= min_area_overlap)
    
    # Download flowlines and VAA for COMIDs in nhd_cats_intrsct
    nhd_flines_intrsct <- nhdplusTools::get_nhdplus(comid = c(nhd_cats_intrsct$featureid),
                                                    realization = 'flowline') %>%
      # reformat variable names to uppercase
      rename_with(.,toupper,id:enabled)
    
    # If omit_zero_area_flines is TRUE, edit df to only return COMIDs where 
    # attribute AREASQKM != 0 (unless there is only one COMID that drains to the 
    # PRMS line, in which case keep the AREASQKM = 0 COMID to make sure that
    # all PRMS/NHM segments are represented in the output data frame).
    if(omit_zero_area_flines & length(nhd_flines_intrsct$COMID) > 1){
      
      nhd_flines_out <- nhd_flines_intrsct %>%
        filter(AREASQKM > 0) 
      
      # Subset COMIDs that would be excluded because AREASQKM is zero, and print
      # a message to the console letting the user know if any of the excluded COMIDs
      # have a reach length greater than zero_area_flines_max_length.
      nhd_flines_excluded <- setdiff(nhd_flines_intrsct,nhd_flines_out)
      
      if(any(nhd_flines_excluded$LENGTHKM > zero_area_flines_max_length)){
        message(sprintf(
          "The following NHDPlusv2 COMIDs are longer than %s km but are being omitted because AREASQKM == 0: \n\n%s\n",
          zero_area_flines_max_length,
          paste(nhd_flines_excluded$COMID[nhd_flines_excluded$LENGTHKM > zero_area_flines_max_length], collapse="\n")))
      }
      
    } else {
      nhd_flines_out <- nhd_flines_intrsct
    }
    
    # Save intersecting catchment COMID's
    comids_all <- data.frame(PRMS_segid = prms_line$subsegid,
                             comid_cat = paste(sort(unique(nhd_flines_out$COMID)),collapse=";"))
    
  } else {
    comids_all <- data.frame(PRMS_segid = prms_line$subsegid,
                             comid_cat = NA)
  }
  
  return(comids_all)
  
}



get_upstream_nhd_catchments <- function(nhd_lines, xwalk_table, prms_line,
                                        omit_zero_area_flines,
                                        zero_area_flines_max_length){
  #' @description This function finds the NHDPlusV2 catchments that directly drain to 
  #' a given segment by navigating the NHDPlusV2 value-added attributes tables using 
  #' the information in the NHDPlus/PRMS reach-to-segment xwalk table
  #' 
  #' @param prms_line sf linestring containing a target PRMS segment
  #' @param xwalk_table data frame containing the NHDPlusV2/PRMS reach-to-segment 
  #' crosswalk information
  #' @param nhd_lines sf object containing NHDPlusV2 flowlines for area of interest
  #' nhd_lines must contain variables COMID, PATHLENGTH, LENGTHKM, and HYDROSEQ
  #' @param omit_zero_area_flines logical; if TRUE, only return NHDPlusv2 reaches 
  #' where the NHDPlus atttribute AREASQKM is greater than zero. 
  #' @param zero_area_flines_max_length integer; zero-area reaches are expected to be
  #' relatively short, so a message is printed to the console if `omit_zero_area_flines` 
  #' is TRUE and `LENGTHKM` of an omitted reach exceeds `zero_area_flines_max_length`.
  #' 
    
  # subset xwalk for prms_line
  xwalk_prms_line <- xwalk_table %>%
    filter(PRMS_segid == prms_line$subsegid)
  
  # Find all tributaries upstream of comid_down
  nhd_reach_down_UT <- nhdplusTools::get_UT(nhd_lines,xwalk_prms_line$comid_down)
  
  # Find Geospatial Fabric (GFv1) POI's (=respective comid_down's) that are upstream 
  # of the PRMS segment of interest
  ntw_POI <- xwalk_table %>%
    filter(PRMS_segid != prms_line$subsegid) %>%
    filter(comid_down %in% nhd_reach_down_UT) %>%
    select(-comid_seg) 
  
  # Find all upstream tributaries for each upstream GF POI
  ntw_POI_UT <- lapply(ntw_POI$comid_down,nhdplusTools::get_UT,network=nhd_lines) %>%
    do.call("c",.)
  
  # By difference, find the tributaries that contribute directly to prms_line
  nhd_reach_diff <- setdiff(nhd_reach_down_UT,ntw_POI_UT)
  
  # Retrieve NHDPlus attributes for COMIDS that contribute directly to prms_line
  nhd_reach_diff_flines <- nhd_lines %>%
    filter(COMID %in% nhd_reach_diff)
  
  # If omit_zero_area_flines is TRUE, edit df to only return COMIDs where 
  # attribute AREASQKM != 0 (unless there is only one COMID that drains to the 
  # PRMS line, in which case keep the AREASQKM = 0 COMID to make sure that
  # all PRMS/NHM segments are represented in the output data frame).
  if(omit_zero_area_flines & length(nhd_reach_diff_flines$COMID) > 1){
    
    nhd_reach_diff_flines_out <- nhd_reach_diff_flines %>%
      filter(AREASQKM > 0) 
    
    # Subset COMIDs that would be excluded because AREASQKM is zero, and print
    # a message to the console letting the user know if any of the excluded COMIDs
    # have a reach length greater than zero_area_flines_max_length.
    nhd_reach_diff_flines_excluded <- setdiff(nhd_reach_diff_flines,nhd_reach_diff_flines_out)
    
    if(any(nhd_reach_diff_flines_excluded$LENGTHKM > zero_area_flines_max_length)){
      message(sprintf(
        "The following NHDPlusv2 COMIDs are longer than %s km but are being omitted because AREASQKM == 0: \n\n%s\n",
        zero_area_flines_max_length,
        paste(nhd_reach_diff_flines_excluded$COMID[nhd_reach_diff_flines_excluded$LENGTHKM > zero_area_flines_max_length], collapse="\n")))
    }
      
  } else {
    nhd_reach_diff_flines_out <- nhd_reach_diff_flines 
  }
  
  # Save contributing catchment COMID's
  comids_all <- data.frame(PRMS_segid = prms_line$subsegid,
                           comid_cat = paste(sort(nhd_reach_diff_flines_out$COMID),collapse=";"))
  
  return(comids_all)
  
}
