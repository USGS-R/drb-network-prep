#' Calculate distance matrix for reaches in network
#'
#' @description Uses igraph to compute a matrix of distances in meters.
#' 
#' @details This function originally developed in the delaware-model-prep repo
#' and modified here for use with NHDPlusv2.
#' https://github.com/USGS-R/delaware-model-prep/blob/main/1_network/src/calc_distance_functions.R
#'
#' @param network sf linestring object representing the network from which to
#' calculate distances. Must contain columns comid, tonode, and fromnode at minimum.
#' 
#' @return distance matrix
#' 
calc_dist_matrices_nhd <- function(network){

  # Prep network (must contain columns comid, tonode,fromnode)
  network <- nhdplusTools::get_tocomid(network, add = TRUE) %>%
    rename_with(.,toupper,comid:enabled) %>%
    select(ID = COMID, toID = TOCOMID, LENGTHKM) %>%
    mutate(length_m = LENGTHKM*1000)
    
  # Code below modified from delaware-model-prep: 1_network/src/calc_distance_functions.R
  # combine reach information so that the "edges" of the river network (the river reaches)
  # become the vertices of the igraph object, where the distance from one vertex to the 
  # next is the distance along the second vertex (reach).
  edges <- network %>%
    st_drop_geometry() %>%
    select(ID, toID, length_m) %>%
    rename(from_reach=ID, to_reach=toID, reach_length = length_m) %>%
    # remove the rows where a COMID empties to 0 (the river outlets); the COMIDs are still present in the network as toID's
    filter(to_reach != 0)
  
  # df must contain a "symbolic edge list in the first two columns.
  # Additional columns are considered as edge attributes"
  graph <-  igraph::graph_from_data_frame(edges, directed = TRUE) 
  
  # calculate symmetric distance:
  dists_complete <- igraph::distances(graph, weights = igraph::edge.attributes(graph)$reach_length, mode='all')
  # Calculate shortest paths FROM each vertex:
  dists_downstream <- igraph::distances(graph, weights = igraph::edge.attributes(graph)$reach_length, mode='out')
  # Calculate shortest paths TO each vertex (flowing upstream only):
  dists_upstream <- igraph::distances(graph, weights = igraph::edge.attributes(graph)$reach_length, mode='in')
  dists_updown <- dists_downstream
  for(i in 1:nrow(dists_downstream)) {
    for(j in 1:ncol(dists_downstream)) {
      if(is.infinite(dists_downstream[i,j]) & !is.infinite(dists_upstream[i,j])) {
        dists_updown[i,j] <- -dists_upstream[i,j]
      }
    }
  }
  out <- list(complete = dists_complete, 
              downstream = dists_downstream, 
              upstream = dists_upstream, 
              updown = dists_updown)
  
  return(out)
}
