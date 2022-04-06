library(targets)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse","sf","purrr", 
                            "sbtools","nhdplusTools")) 

source("1_fetch.R")
source("2_process.R")

dir.create("1_fetch/out/", showWarnings = FALSE)
dir.create("1_fetch/log/", showWarnings = FALSE)
dir.create("2_process/out/", showWarnings = FALSE)
dir.create("2_process/log/", showWarnings = FALSE)

# Define dataset of interest for the national geospatial fabric (used to fetch PRMS catchment polygons):
gf_data_select <- 'GeospatialFabric_National.gdb.zip'

# Define which columns related to HRU polygons to keep
# In national GF, hru_id is identical to hru_id_reg 
gf_cols_select <- c("hru_id_nat","hru_id","hru_segment","region","Shape_Area")

# Define minor HUCs (hydrologic unit codes) that make up the DRB to use in calls to dataRetrieval functions
# Lower Delaware: 0204 subregion (for now, exclude New Jersey Coastal (https://water.usgs.gov/GIS/huc_name.html)
drb_huc8s <- c("02040101","02040102","02040103","02040104","02040105","02040106",
               "02040201","02040202","02040203","02040204","02040205","02040206","02040207")

# Define PRMS segments that require special handling for the Geospatial Fabric (GFv1) to NHDplusV2 crosswalk target
# These segments were flagged because the cumulative NHD catchment areas are considerably different from the 
# summed area of the corresponding HRU's from the NHGFv1. For these segments, contributing NHDPlusV2 catchments 
# are identified using a spatial join with the corresponding PRMS HRU's. More details can be found in the 
# USGS-R/drb-inland-salinity-ml repo: https://github.com/USGS-R/drb-inland-salinity-ml/issues/45
drb_segs_spatial <- c("31_1","135_1","233_1","236_1","244_1","249_1","332_1","354_1","355_1","358_1","396_1","593_1","1256_1","2137_1",
                      "2138_1","2757_1","2758_1","2759_1","2765_1","2766_1","2767_1","2772_1","2797_1")

# Define the PRMS segments that were split within the USGS-R/delaware-model-prep pipeline; these 
# segments require special handling to return proper catchment areas in place of NHGFv1 HRUs
GFv1_segs_split <- c("3_1","3_2","8_1","8_2","51_1","51_2")


# Return the complete list of targets
c(p1_targets_list, p2_targets_list)


