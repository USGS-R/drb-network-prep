Active development moved to https://code.usgs.gov/wma/wp/drb-network-prep

# drb-network-prep
Code repo to prepare the river network for modeling in the Delaware River Basin  
  
## Prep the DRB river network for modeling  
Nothing here yet. Functions to prep the DRB network are currently included in the delaware-model-prep [repo](https://github.com/USGS-R/delaware-model-prep).  
  
## Create NHGF-NHDv2 crosswalk table
Functions here create a crosswalk table between PRMS segments (currently, NHGFv1) and NHDPlusV2 flowlines.  
  
The xwalk table was initially created as part of the DRB inland salinity project, and the full commit history can be referenced within that project's [repo](https://github.com/USGS-R/drb-inland-salinity-ml). The xwalk table contains the `PRMS_segid` as well as columns `comid_down`, which represents the NHDPlusV2 COMID at the downstream end of each PRMS segment, `comid_seg`, which represents the NHDPlusV2 COMID's that intersect each PRMS segment, and `comid_cat`, which represents all of the NHDPlusV2 COMID's that directly drain to the PRMS segment, including headwater tributaries. For reference, this repo also includes two log files that contain metadata for each target built in the pipeline as well as the sf external dependency versions that were used to build the xwalk table.  


