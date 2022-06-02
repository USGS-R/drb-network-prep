create_GFv1_NHDv2_xwalk <- function(prms_lines, nhd_lines, prms_hrus,
                                    min_area_overlap, drb_segs_spatial,
                                    omit_divergences = FALSE, 
                                    omit_zero_area_flines = FALSE){
  #' @description This function outputs a table to facilitate cross-walking between
  #'  PRMS river segments (from the Geospatial Fabric, GFv1) and NHDPlusV2 flowlines 
  #'  and catchments for the Delaware River Basin (DRB). 
  #' 
  #' @param prms_lines sf object containing the PRMS river segments for the DRB
  #' @param nhd_lines sf object containing NHDPlusV2 flowlines for area of interest
  #' nhd_lines must contain variables COMID, PATHLENGTH, LENGTHKM, HYDROSEQ, STREAMORDE, 
  #' STREAMCALC, FROMNODE,and TONODE
  #' @param prms_hrus sf (multi)polygon containing the HRU's from the GFv1
  #' @param min_area_overlap float; value indicating the minimum proportion of NHDPlusV2 
  #' catchment area that overlaps the PRMS polygon in order to be retained
  #' @param drb_segs_spatial character vector containing the identity of PRMS segments 
  #' that require special handling. For these segments, the contributing NHD catchments 
  #' will be determined using a spatial join with the PRMS HRU polygons.
  #' @param omit_divergences logical; if TRUE, only return NHDv2 reaches that are part
  #' of the dendritic river network (i.e., streamcalc = streamorder). Defaults to FALSE.
  #' @param omit_zero_area_flines logical; if TRUE, only return NHDv2 reaches where the 
  #' NHD attribute "AREASQKM" is greater than zero. Defaults to FALSE.
  #
  
  # Pass NHDPlusV2 flowlines with or without divergences to downstream crosswalk functions
  if(omit_divergences){
    nhd_lines_select <- nhd_lines %>%
      filter(STREAMORDE == STREAMCALC)
  } else {
    nhd_lines_select <- nhd_lines
  }
  
  # Add special handling for split segments so that all split segments have a value
  # for segidnat (instead of NA)
  prms_lines_complete <- prms_lines %>%
    group_by(subsegseg) %>%
    mutate(segidnat = segidnat[!is.na(segidnat)]) %>%
    ungroup()
  
  # NHDPlusV2 reaches where AREASQKM equals zero are expected to be relatively short,
  # so assign a max reach length (km) such that a message is printed to the console if 
  # omit_zero_area_flines is TRUE and an omitted reach is longer than the 25 percent
  # quantile of all reach lengths
  zero_area_flines_max_length <- as.numeric(
    quantile(nhd_lines_select$LENGTHKM, 0.25, na.rm = TRUE))

  # find NHDPlusV2 COMID's that intersect PRMS segments
  reach_to_seg_xwalk <- prms_lines_complete %>%
    split(.,.$subsegid) %>%
    purrr::map(.,pair_nhd_reaches,
               nhd_lines = nhd_lines_select, 
               omit_divergences = omit_divergences, 
               omit_zero_area_flines = omit_zero_area_flines,
               zero_area_flines_max_length = zero_area_flines_max_length
               ) %>%
    purrr::map(.,summarize_paired_comids) %>%
    bind_rows()
  
  # find NHDPlusV2 COMID's that drain directly to PRMS segments
  cats_to_seg_xwalk <- prms_lines_complete %>%
    split(.,.$subsegid) %>%
    purrr::map(.,pair_nhd_catchments,
               prms_hrus = prms_hrus,
               min_area_overlap = min_area_overlap,
               nhd_lines = nhd_lines_select,
               xwalk_table = reach_to_seg_xwalk,
               drb_segs_spatial = drb_segs_spatial,
               omit_zero_area_flines = omit_zero_area_flines,
               zero_area_flines_max_length = zero_area_flines_max_length) %>%
    bind_rows() %>%
    left_join(reach_to_seg_xwalk,.,by="PRMS_segid")
  
  return(cats_to_seg_xwalk)
}

