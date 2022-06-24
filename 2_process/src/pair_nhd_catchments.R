#' @description This function finds the NHDPlusv2 catchments that directly drain to 
#' a given segment by using the information in the NHM/NHDPlus reach-to-segment 
#' crosswalk table and navigating the NHDPlusv2 value-added attributes tables 
#' 
#' @param nhd_lines sf object containing NHDPlusV2 flowlines for area of interest
#' nhd_lines must contain variables COMID, PATHLENGTH, LENGTHKM, and HYDROSEQ
#' @param xwalk_table data frame containing the NHM/NHDPlus reach-to-segment 
#' crosswalk information
#' @param prms_line sf linestring containing a target PRMS segment
#' @param omit_zero_area_flines logical; if TRUE, only return NHDPlusv2 reaches 
#' where the NHDPlus atttribute AREASQKM is greater than zero. 
#' @param zero_area_flines_max_length integer; zero-area reaches are expected to be
#' relatively short, so a message is printed to the console if `omit_zero_area_flines` 
#' is TRUE and `LENGTHKM` of an omitted reach exceeds `zero_area_flines_max_length`.
#' 
pair_nhd_catchments <- function(nhd_lines, xwalk_table, prms_line,
                                omit_zero_area_flines, zero_area_flines_max_length){

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
