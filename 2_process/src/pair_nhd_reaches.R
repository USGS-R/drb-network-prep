pair_nhd_reaches <- function(nhd_lines, prms_line, omit_divergences, omit_zero_area_flines){
  #' 
  #' @description Function to pair PRMS segments with associated NHDPlusV2 flowlines
  #'
  #' @param nhd_lines sf object containing NHDPlusV2 flowlines for area of interest
  #' nhd_lines must contain variables COMID,STREAMORDE,STREAMCALC,HYDROSEQ,FROMNODE,TONODE
  #' @param prms_line sf linestring representing the target PRMS segment
  #' @param omit_divergences logical; if TRUE, only return paired NHDv2 reaches that are 
  #' part of the dendritic river network (i.e., streamcalc = streamorder). Defaults to TRUE.
  #' @param omit_zero_area_flines logical; if TRUE, only return paired NHDPlusv2 reaches where
  #' the NHDPlus attribute "AREASQKM" is greater than zero. 
  #'
  #' @value Data frame containing the paired NHDPlusV2 reaches
  
  # Project PRMS and NHD reach lines:
  prms_line_proj <- sf::st_transform(prms_line,5070)
  nhd_lines_proj <- sf::st_transform(nhd_lines,5070)
  
  # Create small buffer (0.3 m) around PRMS segment:
  prms_line_buffer <- sf::st_buffer(prms_line_proj,dist = 0.3)
  
  # Find intersection between NHD reaches and buffered PRMS polygon, and calculate lengths of intersecting NHD reaches
  lines_int <- suppressWarnings(sf::st_intersection(nhd_lines_proj,prms_line_buffer)) %>%
    mutate(len = as.numeric(sf::st_length(.))) %>%
    sf::st_drop_geometry()

  # Find associated NHD reaches that overlap by >5 m (indicating that the NHD and PRMS lines likely overlap rather than just touch)
  nhd_paired <- lines_int %>% 
    filter(len > 5)
  
  # Identify furthest upstream and downstream COMIDs within the group of associated NHD reaches (omitting divergences)
  nhd_paired_down <- nhd_paired %>%
    filter(STREAMORDE==STREAMCALC) %>%
    filter(HYDROSEQ==min(HYDROSEQ))
  
  nhd_paired_up <- nhd_paired %>%
    filter(STREAMORDE==STREAMCALC) %>%
    filter(HYDROSEQ==max(HYDROSEQ))
  
  # Add in special handling for select PRMS segments 
  # PRMS line "341_1" doesn't completely overlap NHD, leading to misleading upstream HYDROSEQ id
  if(prms_line$subsegid=="341_1"){
    nhd_paired_up <- nhd_paired %>%
      filter(COMID == "4480903")
  }
  # The upstream COMID associated with PRMS subsegid "2765_1" likely does not drain to the 
  # downstream COMID, leading to infinite recursion errors if omit_divergences is TRUE
  if(omit_divergences & prms_line$subsegid == "2765_1"){
    nhd_paired_up <- nhd_paired %>%
      filter(COMID == "9484458")
  }
  
  # To make sure we're not missing any NHD reaches between downstream/upstream COMIDs,
  # recursively traverse NHD from upstream COMID until we find downstream COMID
  if(nhd_paired_down$COMID==nhd_paired_up$COMID){
    between_lines <- nhd_paired_down
  }else{
    between_lines <- traverse_nhd(nhd_lines = nhd_lines,
                                  paired_flines = nhd_paired,
                                  down_comid = nhd_paired_down$COMID,
                                  up_comid=nhd_paired_up$COMID)
  }
  
  # Save df containing paired NHD reaches
  df_paired <- between_lines %>%
    mutate(PRMS_segid = prms_line$subsegid) %>%
    # save segidnat column from PRMS lines, adding special handling
    # to impute segidnat values for split segments
    mutate(segidnat = case_when(
      PRMS_segid == "3_1" ~ as.integer("1437"),
      PRMS_segid == "8_1" ~ as.integer("1442"),
      PRMS_segid == "51_1" ~ as.integer("1485"),
      TRUE ~ prms_line$segidnat)) 
  
  # If omit_zero_area_flines is TRUE, edit df to only return COMIDs where 
  # attribute AREASQKM != 0 (unless there is only one COMID matched to the 
  # PRMS line, in which case keep the AREASQKM = 0 COMID to make sure that
  # all PRMS/NHM segments are represented in the output data frame).
  if(omit_zero_area_flines & length(df_paired$COMID) > 1){
    df_out <- df_paired %>%
      filter(AREASQKM > 0) %>%
      select(PRMS_segid,segidnat,COMID,HYDROSEQ,LEVELPATHI,REACHCODE,STREAMORDE,STREAMCALC)
  } else {
    df_out <- df_paired %>%
      select(PRMS_segid,segidnat,COMID,HYDROSEQ,LEVELPATHI,REACHCODE,STREAMORDE,STREAMCALC)
  }
  
  return(df_out)
  
}



traverse_nhd <- function(nhd_lines,paired_flines,down_comid,up_comid){
  #' 
  #' @description Function to traverse NHDPlusV2 between user-specified upstream/downstream COMIDs using FromNode and ToNode attributes
  #'
  #' @param nhd_lines sf object containing NHDPlusV2 flowlines for area of interest
  #' nhd_lines must contain variables COMID, FROMNODE, and TONODE
  #' @param paired_flines sf object containing NHDPlusV2 reaches associated with a target PRMS segment
  #' @param down_comid integer containing the COMID of the most downstream NHDPlusV2 reach associated with a target PRMS segment
  #' @param up_comid integer containing the COMID of the most upstream NHDPlusV2 reach associated with a target PRMS segment
  #'
  #' @value data frame containing all of the NHDPlusV2 reaches between down_comid and up_comid
  
  nhd_lines_df <- nhd_lines %>%
    sf::st_drop_geometry() 
  
  # Identify starting (most upstream) ToNode
  start_up <- nhd_lines_df %>% 
    filter(COMID %in% up_comid)
  up_tonode <- start_up$TONODE
  
  # Recursively find downstream ToNode
  down <- nhd_lines_df %>% 
    filter(FROMNODE %in% up_tonode)
  
  # If there are multiple downstream segments (i.e., at a confluence), give it a hint and filter for the COMID within paired_flines
  if(length(down$COMID) > 1){
    down <- down %>% filter(COMID %in% paired_flines$COMID)
  }
  
  if(down_comid %in% down$COMID){
    net <- rbind(start_up,down)
    return(net)
  } else {
    net <- rbind(start_up,down,traverse_nhd(nhd_lines,paired_flines,down_comid,down$COMID)) %>% 
      filter(!duplicated(COMID))
  }
  
}



summarize_paired_comids <- function(paired_nhd_df){
  #' 
  #' @description Function to summarize the paired NHDPlusV2 reaches for each PRMS segment 
  #'
  #' @param paired_nhd_df data frame containing the paired NHDPlusV2 reaches associated with each PRMS segment,
  #' where each row represents one PRMS segment. comid_down contains the NHDPlusV2 reach at the downstream
  #' end of the PRMS segment, whereas comid_seg represents all of the contributing NHDPlusV2 reaches.
  
  comids_out <- paired_nhd_df %>%
    group_by(PRMS_segid) %>% 
    # concatenate all paired NHDPlusV2 reaches into one column, comid_seg
    mutate(comid_seg = paste(sort(unique(COMID)),collapse=";")) %>%
    ungroup() %>%
    # identify most downstream NHDPlusV2 reach (comid_down) for each PRMS segment that does not represent a divergence
    filter(STREAMORDE==STREAMCALC) %>%
    filter(HYDROSEQ==min(HYDROSEQ)) %>% 
    select(PRMS_segid,segidnat,COMID,comid_seg) %>%
    rename(comid_down = COMID)
  
  return(comids_out)
  
}
