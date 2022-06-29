library(targets)

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c("tidyverse","sf", "purrr", "sbtools",
                            "nhdplusTools", "igraph", "reticulate")) 

source("1_fetch.R")
source("2_process.R")

dir.create("1_fetch/out/", showWarnings = FALSE)
dir.create("1_fetch/log/", showWarnings = FALSE)
dir.create("2_process/out/", showWarnings = FALSE)
dir.create("2_process/log/", showWarnings = FALSE)

# Define dataset of interest from the national hydrologic geospatial fabric, NHGF 
# (used to fetch PRMS catchment polygons); values entered below can take two forms 
# to indicate preference for the national-scale data ("GeospatialFabric_National.gdb.zip) 
# or the regional-scale data ("GeospatialFabricFeatures_XX.zip"), where "XX" refers to 
# the region id/vpu of interest.
gf_data_select <- 'GeospatialFabric_National.gdb.zip'

# Define which columns related to HRU polygons to keep
# In national GF, hru_id is identical to hru_id_reg 
gf_cols_select <- c("hru_id_nat","hru_id","hru_segment","region","Shape_Area")

# Define minor HUCs (hydrologic unit codes) that make up the DRB to use in calls to dataRetrieval functions
# Lower Delaware: 0204 subregion (for now, exclude New Jersey Coastal (https://water.usgs.gov/GIS/huc_name.html)
drb_huc8s <- c("02040101","02040102","02040103","02040104","02040105","02040106",
               "02040201","02040202","02040203","02040204","02040205","02040206","02040207")

# Define the PRMS segments that were split within the USGS-R/delaware-model-prep pipeline; these 
# segments require special handling to return proper catchment areas in place of NHGFv1 HRUs
GFv1_segs_split <- c("3_1","3_2","8_1","8_2","51_1","51_2")


# Return the complete list of targets
c(p1_targets_list, p2_targets_list)


