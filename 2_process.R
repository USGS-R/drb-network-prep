source("2_process/src/pair_nhd_reaches.R")
source("2_process/src/pair_nhd_catchments.R")
source("2_process/src/create_GFv1_NHDv2_xwalk.R")
source("2_process/src/munge_GFv1_catchments.R")
source("2_process/src/write_data.R")
source("2_process/src/write_ind_files.R")

p2_targets_list <- list(
  
  # Pair PRMS segments (NHGFv1) with intersecting NHDPlusV2 reaches and contributing 
  # NHDPlusV2 catchments
  # 1) Crosswalk table based on full NHD network
  tar_target(
    p2_prms_nhdv2_xwalk,
    create_GFv1_NHDv2_xwalk(prms_lines = p1_GFv1_reaches_sf,
                            nhd_lines = p1_nhdv2reaches_sf,
                            # if TRUE, only return NHDv2 reaches that are part of the
                            # dendritic river network
                            omit_divergences = FALSE, 
                            # omit NHDv2 reaches where attribute AREASQKM is equal to 0; 
                            # these reaches may represent small extensions, sometimes
                            # artificial paths, or flowlines through a waterbody
                            omit_zero_area_flines = FALSE)
  ),
  # 2) Crosswalk table with divergent reaches omitted (i.e., return paired COMIDs that 
  # belong to the dendritic river network)
  tar_target(
    p2_prms_nhdv2_xwalk_omit_divergences,
    create_GFv1_NHDv2_xwalk(prms_lines = p1_GFv1_reaches_sf,
                            nhd_lines = p1_nhdv2reaches_sf,
                            omit_divergences = TRUE, 
                            omit_zero_area_flines = FALSE)
  ),
  # 3) Crosswalk table with zero-area flowlines removed (i.e., don't return any COMIDs
  # where the AREASQKM attribute equals zero). Divergences are retained to avoid any
  # "holes" in the catchment areas. See https://github.com/USGS-R/drb-network-prep/issues/18
  tar_target(
    p2_prms_nhdv2_xwalk_omit_zero_area,
    create_GFv1_NHDv2_xwalk(prms_lines = p1_GFv1_reaches_sf,
                            nhd_lines = p1_nhdv2reaches_sf,
                            omit_divergences = FALSE, 
                            omit_zero_area_flines = TRUE)
  ),
  
  # Save GFv1-NHDv2 xwalk table 
  # 1) Crosswalk table based on full NHD network
  tar_target(
    p2_prms_nhdv2_xwalk_csv,
    write_to_csv(p2_prms_nhdv2_xwalk,"2_process/out/GFv1_NHDv2_xwalk.csv"),
    format = "file"
  ),
  # 2) Crosswalk table with divergent reaches omitted
  tar_target(
    p2_prms_nhdv2_xwalk_omit_divergences_csv,
    write_to_csv(p2_prms_nhdv2_xwalk_omit_divergences,"2_process/out/GFv1_NHDv2_xwalk_dendritic.csv"),
    format = "file"
  ),
  # 3) Crosswalk table with zero-area flowlines removed
  tar_target(
    p2_prms_nhdv2_xwalk_omit_zero_area_csv,
    write_to_csv(p2_prms_nhdv2_xwalk_omit_zero_area,"2_process/out/GFv1_NHDv2_xwalk_omit_zero_area.csv"),
    format = "file"
  ),
  
  # Reshape full GFv1-NHDv2 xwalk table to return all COMIDs that drain to each PRMS segment
  tar_target(
    p2_drb_comids_all_tribs, 
    p2_prms_nhdv2_xwalk %>%
      select(PRMS_segid, comid_cat) %>% 
      tidyr::separate_rows(comid_cat,sep=";") %>% 
      rename(COMID = comid_cat)
  ),
  
  # Process catchments so that every PRMS segment has >= 1 corresponding HRU; In addition,
  # adjust catchments for 3 segments that were split in delaware-model-prep pipeline
  # https://github.com/USGS-R/delaware-model-prep
  tar_target(
    p2_GFv1_catchments_edited_sf,
    munge_GFv1_catchments(prms_lines = p1_GFv1_reaches_sf,
                          prms_hrus = p1_GFv1_catchments_sf,
                          segs_w_comids = p2_drb_comids_all_tribs,
                          segs_split = GFv1_segs_split,
                          crs_out = 5070,
                          verbose = TRUE)
  ),
  
  # Save processed GFv1 catchments as a geopackage
  tar_target(
    p2_GFv1_catchments_edited_gpkg,
    write_sf(p2_GFv1_catchments_edited_sf,
             dsn = "2_process/out/GFv1_catchments_edited.gpkg", 
             layer = "GFv1_catchments_edited", 
             driver = "gpkg",
             quiet = TRUE,
             # overwrite layer if already exists
             append = FALSE),
    format = "file"
  ),
  
  # Create and save indicator file
  # argument force_dep must contain name of an upstream target to force dependencies
  # and build this target when a log file already exists
  tar_target(
    p2_data_summary_csv,
    write_ind_files("2_process/log/GFv1_data_summary.csv",
                    force_dep = c(p2_prms_nhdv2_xwalk_csv,
                                  p2_prms_nhdv2_xwalk_omit_divergences_csv,
                                  p2_prms_nhdv2_xwalk_omit_zero_area_csv),
                    target_names = c("p1_GFv1_reaches_sf","p1_GFv1_catchments_sf","p1_nhdv2reaches_sf",
                                     "p1_nhdv2_catchments_sf","p1_nhdv2_catchments_gpkg",
                                     "p2_prms_nhdv2_xwalk","p2_prms_nhdv2_xwalk_csv",
                                     "p2_prms_nhdv2_xwalk_omit_divergences",
                                     "p2_prms_nhdv2_xwalk_omit_divergences_csv",
                                     "p2_prms_nhdv2_xwalk_omit_zero_area",
                                     "p2_prms_nhdv2_xwalk_omit_zero_area_csv")),
    format = "file"),
  
  # Create and save sf session info
  tar_target(
    p2_sf_version_csv,
    write_session_info("2_process/log/sf_version_info.csv"),
    format = "file"
  )
  
)






