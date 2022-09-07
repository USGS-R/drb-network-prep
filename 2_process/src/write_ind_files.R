write_ind_files <- function(fileout,target_names,force_dep) {
  #' 
  #' @description Function to save indicator files to track data changes over time
  #'
  #' @param fileout a character string that indicates the name of the file to be saved, including path and file extension 
  #' @param target_names a character string or vector of strings containing the target names of interest
  #' @param force_dep name of an upstream target; force_dep does not do anything in the function and so
  #' is only acting as a placeholder to force dependencies when building the metadata table
  #'
  #' @value Returns a csv file containing the target metadata
  
  # Create indicator table
  ind_tbl <- tar_meta(all_of(target_names)) %>%
    select(tar_name = name, 
           latest_build_date = time, 
           filepath = path, 
           hash = data) %>%
    mutate(filepath = unlist(filepath),
           latest_build_date = as.Date(latest_build_date))
  
  # Save indicator table
  readr::write_csv(ind_tbl, fileout)
  return(fileout)
}


write_session_info <- function(fileout){
  
  # save R version info and select package versions
  version_info <- data.frame(Soft = c("R","sf","nhdplusTools","Contact:"),
                             version = c(R.version$version.string, 
                                         as.character(packageVersion("sf")), 
                                         as.character(packageVersion("nhdplusTools")),
                                         Sys.info()["user"]))

  # save sf external dependency versions
  extSoft_version_info <- sf::sf_extSoftVersion() %>% 
    data.frame(Soft = names(.),version = .)
  
  # combine version info
  version_info_combined <- extSoft_version_info %>%
    bind_rows(.,version_info) 
  
  # Save session info table
  write_to_csv(version_info_combined,fileout)
  
  return(fileout)

}

