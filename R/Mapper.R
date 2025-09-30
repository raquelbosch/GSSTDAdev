#' @title one_D_Mapper
#'
#' @description Wrapping function to carry out the complete process.
#'
#' @param mapper_object_ini Mapper TDA initializated object generated
#' by \code{mapper} function.
#'
#' @export
#' @return A \code{mapper_obj} object. It contains the values of the intervals
#' (interval_data), the samples included in each interval (sample_in_level),
#' information about the cluster to which the individuals in each interval
#' belong (clustering_all_levels), a list including the individuals contained
#' in each detected node (node_samples), a list including the nodes contained in
#' each interval (nodes_in_lev), their size (node_sizes), the average of the
#' filter function values of the individuals of each node (node_average_filt)
#' and the adjacency matrix linking the nodes (adj_matrix).
#' Moreover, information is provided on the number of nodes, the average node
#' size, the standard deviation of the node size, the number of connections
#' between nodes, the proportion of connections to all possible connections
#' and the number of ramifications.
one_D_Mapper <- function(mapper_object_ini){

  data <- mapper_object_ini[["data"]]
  filter_values <- mapper_object_ini[["filter_values"]]

  #Getting intervals.
  if (mapper_object_ini[["type_covering"]] == "uniform"){
    interval_data <- get_intervals_One_D_uniform_covering(filter_values,
                                                          mapper_object_ini[["num_intervals"]],
                                                          mapper_object_ini[["percent_overlap"]])
  } else {
    interval_data <- get_intervals_One_D(filter_values, mapper_object_ini[["num_intervals"]],
                                         mapper_object_ini[["percent_overlap"]])
  }

  #Getting samples on each interval.
  samp_in_lev <- samples_in_levels(interval_data, filter_values)

  #Clustering all levels.
  test_clust_all_levels <- clust_all_levels(data, samp_in_lev,
                                            mapper_object_ini[["distance_type"]],
                                            mapper_object_ini[["clustering_type"]],
                                            mapper_object_ini[["linkage_type"]],
                                            mapper_object_ini[["optimal_clustering_mode"]],
                                            mapper_object_ini[["silhouette_threshold"]],
                                            mapper_object_ini[["num_bins_when_clustering"]],
                                            mapper_object_ini[["dim_reduction"]])
  #Transforming levels into nodes.
  node_samples <- levels_to_nodes(test_clust_all_levels)

  # Getting nodes on each interval.
  nodes_in_lev <- nodes_in_levels(test_clust_all_levels)

  #Computing adjacency matrix.
  adj_matrix_out <- compute_node_adjacency(node_samples)

  node_sizes = unlist(lapply(node_samples,length))
  # average of the filter function of each node
  node_average_filt = lapply(node_samples,function(x,y) mean(y[x]),filter_values)

  # additional parameters
  n_nodes <- length(node_sizes)
  av_node_size <- mean(node_sizes)
  sd_node_size <- stats::sd(node_sizes)

  adj_mat <- adj_matrix_out
  upper_tri <- upper.tri(adj_mat)
  lower_tri <- lower.tri(adj_mat)

  n_connections <- sum(adj_mat[upper_tri] == 1)
  prop_connections <- n_connections/length(adj_mat[upper_tri])

  adj_mat[lower_tri] <- t(adj_mat)[lower_tri]
  diag(adj_mat) <- 0
  n_ramifications <- colSums(adj_mat)-2
  n_ramifications[n_ramifications %in% c(-1,-2)] <- 0
  n_ramifications <- sum(n_ramifications)

  #Generating the object of the output data
  mapper_object <- list("interval_data" = interval_data,
                        "sample_in_level" = samp_in_lev,
                        "clustering_all_levels" = test_clust_all_levels,
                        "node_samples" = node_samples,
                        "nodes_in_lev" = nodes_in_lev,
                        "node_sizes" = node_sizes,
                        "node_average_filt" = node_average_filt,
                        "adj_matrix" = adj_matrix_out,
                        "n_sizes" = n_nodes,
                        "average_nodes"= av_node_size,
                        "standard_desviation_nodes" = sd_node_size,
                        "number_connections" = n_connections,
                        "proportion_connections" = prop_connections,
                        "number_ramifications" = n_ramifications)

  class(mapper_object) <- "mapper_object"
  return(mapper_object)
}

