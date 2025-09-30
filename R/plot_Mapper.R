#' @title Map to color
#' @description
#' Auxiliary function that maps a numeric vector, the average node
#' filtering function values, to a color vector.
#' @param x A vector of numeric values storing the average filtering
#' function values found in the samples placed into a specific node.
#' @param limits A two element numeric vector including the range of values.
#' This is optional.
#' @return A vector of the same length of x with colors ranging from blue to
#' red.
#' @import grDevices
map_to_color <- function(x, limits=NULL){
  pallette_ob <-  grDevices::colorRampPalette(colors = c("blue","red"))(100)
  if(is.null(limits)){
    limits=range(x)}
  map_to_col <- pallette_ob[base::findInterval(x,base::seq(limits[1],limits[2],
                                                           length.out=length(pallette_ob)+1),
                                               all.inside=TRUE)]
  return(map_to_col)
}

#' @title list_to_df
#' @description Auxiliary function. It transforms a list into a two-column
#' data frame. The first column contains the names of the elements in the list
#' and the second column contains the values of each element in the list.
#' Each row of the resulting data.frame contains each element of each vector
#' that was part of the list with the name that said vector received in the
#' list.
#' @param list_obj List of vectors that will be transformed into a data.frame.
#' The elements of the list must be named.
#' @param columns_names Optional. Names of the columns in the output data frame.
#' @return Output data.frame. Each row contains each element of each vector
#' that was part of the list with the name that said vector received in the list.
list_to_df <- function(list_obj, columns_names = c("col1", "col2")){
  df_obj <- do.call(rbind, lapply(names(list_obj), function(name) {
    data.frame(col1 = name, col2 = list_obj[[name]],
               stringsAsFactors = FALSE)
  }))
  colnames(df_obj)<- columns_names
  return(df_obj)
}

#' @title percent_table
#' @description Auxiliary function. It creates a contingency table with the
#' percentages of the categories of two variable of interest (columns) in each
#' of the nodes of the mapper graph (rows).
#' @param var_1 First variable, corresponding to the rows in the contingency
#' table
#' @param var_2 Second variable, corresponding to the columns in the contingency
#' table
#' @return Contingency table of the two variables in percentage.
percent_table <- function(var_1, var_2){
  table_vals <- table(var_1, var_2)
  row_summed <- rowSums(table_vals)
  table_vals <- (table_vals/row_summed)*100
  return(table_vals)
}

#' @title percent_table_nodes
#' @description It creates a contingency table showing the relative frequencies of
#' the variable of interest in each of the nodes of the mapper
#' graph.
#' @param mapper_object A list produced as an output of the \code{one_D_Mapper}
#' function.
#' @param var_of_interest Vector of type "character" or "factor" that collects
#' the variable of interest. It must have as names the names of the samples
#' contained in the mapper object.
#' @return Contingency table in percentages showing the relative frequencies of
#' the variable of interest (columns) in each of the nodes of the mapper
#' graph (rows).
percent_table_nodes <- function(mapper_object, var_of_interest) {
  # Create data.frame
  df_node_sample <- list_to_df(mapper_object$node_samples, c("Node", "Sample"))
  df_node_sample$var_of_interest <- var_of_interest[match(df_node_sample$Sample,
                                                          names(var_of_interest))]
  # Create contingency table
  table_node_var <- percent_table(df_node_sample$Node,
                                  df_node_sample$var_of_interest)
  table_node_var <- table_node_var[paste("Node", 1:length(unique(df_node_sample$Node)),
                                         sep = "_"),]
  return(table_node_var)
}

#' @title Plot mapper
#' @description This function produces an interactive network plot using
#' the \code{visNetork} function from the mapper results. By default, the graph
#' is coloured with the average value of the filtering function at each node.
#' If a value of \code{variable_by_colour} is provided, the nodes will be
#' coloured by the average value of that variable at each node. Nodes with
#' lower average values are coloured blue and those with higher values are
#' coloured red.
#' @param mapper_object A list produced as an output of the \code{one_D_Mapper}
#' function.
#' @param trans_node_size Logical, it indicates whether you want to resize
#' the size of the nodes. \code{TRUE} default option.
#' @param exp_to_res Only necessary if trans_node_size is \code{TRUE}. An
#' exponent of the form 1/n to which the node sizes must be raised in order
#' to resize them.
#' @param variable_by_color Optional. Numeric vector containing the values of
#' the variable by which the graph is to be coloured. It must have the same
#' names as those contained in the mapper_object.
#' @return Plots an interactive network using the \code{visNetwork} function.
#' @export
#' @import visNetwork
#' @examples
#' \donttest{
#' # Create data object
#' data_object <- list("full_data" = full_data, "case_tag" = case_tag)
#' class(data_object) <- "data_object"
#'
#' #Select gene from data object
#' gene_selection_object <- gene_selection(data_object, survival_time,
#'                                         survival_event,
#'                                         gene_select_surv_type = "Top_Bot",
#'                                         percent_gen_select_for_fun_filt = 1,
#'                                         gene_select_mapper_metric = "mad",
#'                                         percent_gen_select_for_mapper = 5)
#'
#' mapper_object <- mapper(data = gene_selection_object[["case_genes_disease_component"]],
#' filter_values = gene_selection_object[["filter_values"]],
#' num_intervals = 5, percent_overlap = 40,
#' type_covering = "uniform", distance_type = "correlation",
#' clustering_type = "hierarchical",
#' linkage_type = "single")
#' plot_mapper(mapper_object)}
plot_mapper <- function(mapper_object, trans_node_size = TRUE, exp_to_res = 1/2,
                        variable_by_color = NULL){

  if (is.null(variable_by_color)) {
    color_vector <- map_to_color(base::unlist(mapper_object[["node_average_filt"]]))
  } else {
    check_plot_vectors(mapper_object, variable_by_color,
                       class_vector = "numeric")
    # Calculate the average:
    var_by_color_average <- lapply(mapper_object$node_samples,
                                   function(x,y) mean(y[x]), variable_by_color)
    color_vector <- map_to_color(var_by_color_average)
  }

  arr_ind <- base::which(arr.ind = TRUE, mapper_object[["adj_matrix"]] == 1)
  df_out <- base::data.frame(base::rownames(mapper_object[["adj_matrix"]])[arr_ind[,1]],
                             base::colnames(mapper_object[["adj_matrix"]])[arr_ind[,2]])
  df_out <- base::cbind(arr_ind, df_out)
  base::rownames(df_out) <- 1:base::nrow(df_out)
  base::colnames(df_out) <- c("from","to","from_name","to_name")
  nodes_to_net <- base::unique(base::data.frame(c(df_out[,1]-1,df_out[,2]-1),
                                                c(df_out[,3],df_out[,4])))
  nodes_to_net$node_size <- mapper_object[["node_sizes"]]
  nodes_to_net$ori_size <- mapper_object[["node_sizes"]] #node size with trans_node_size
  if(trans_node_size){
    nodes_to_net$node_size <- (nodes_to_net$node_size)^exp_to_res
  }
  base::colnames(nodes_to_net) <- c("id", "label", "size", "ori_size")
  nodes_to_net$color <- color_vector

  edges_to_net <- df_out[,c(1,2)]-1
  base::colnames(edges_to_net) <- c("from","to")

  nodes <- data.frame(
    id = nodes_to_net$id,
    label = nodes_to_net$label,
    size = nodes_to_net$size,  # Example sizes
    title = paste("Node Size:", nodes_to_net$ori_size),  # Tooltip content
    color = nodes_to_net$color,
    stringsAsFactors = FALSE  # Ensure string handling
  )

  visNetwork::visNetwork(nodes,edges_to_net[!edges_to_net$from == edges_to_net$to,],) %>%
    visNetwork::visNodes(shape = "dot", scaling = list(label = list(enabled = TRUE, min = 10, max = 30))) %>%
    visNetwork::visEdges(smooth = TRUE) %>%
    visNetwork::visLayout(randomSeed = 2) %>%
    visNetwork::visPhysics(solver = "forceAtlas2Based")
}

#' @title Plot Mapper with Pie Charts
#' @description It reproduces the graph of a Mapper object in which each node
#' is represented by a pie chart showing the proportion of each category of a
#' variable of interest in each node.
#' @param mapper_object A list produced as an output of the \code{one_D_Mapper}
#' function.
#' @param var_of_interest Vector of type "character" or "factor" that collects
#' the variable of interest. It must have as names the names of the samples
#' contained in the mapper object. Therefore, within a GSSTDA, we recommend that
#' you first select only those items corresponding to pathological samples.
#' @param name_var_of_interest Name of the variable of interest to be displayed
#' in the legend. "Variable of interest" by default.
#' @param show_node_names Option to indicate whether node labels should be
#' displayed (\code{TRUE}) or not (\code{FALSE}). By default, it takes the value
#' \code{FALSE}. Any value other than \code{FALSE} will be interpreted as
#' \code{TRUE}.
#' @return Mapper plot using \code{igraph} in which each node
#' is represented by a pie chart showing the proportion of each category of a
#' variable of interest in each node.
#' @export
#' @importFrom igraph graph_from_adjacency_matrix
#' @importFrom igraph V
#' @importFrom graphics legend
#' @import RColorBrewer
#' @examples
#' \donttest{
#' # Create data object
#' data_object <- list("full_data" = full_data, "case_tag" = case_tag)
#' class(data_object) <- "data_object"
#'
#' #Select gene from data object
#' gene_selection_object <- gene_selection(data_object, survival_time,
#'                                         survival_event,
#'                                         gene_select_surv_type = "Top_Bot",
#'                                         percent_gen_select_for_fun_filt = 1,
#'                                         gene_select_mapper_metric = "mad",
#'                                         percent_gen_select_for_mapper = 5)
#'
#' mapper_object <- mapper(data = gene_selection_object[["case_genes_disease_component"]],
#' filter_values = gene_selection_object[["filter_values"]],
#' num_intervals = 5, percent_overlap = 40,
#' type_covering = "uniform", distance_type = "correlation",
#' clustering_type = "hierarchical",
#' linkage_type = "single")
#'
#' # We want to represent the percentage of mortalities for each node.
#' names(survival_event) <- colnames(full_data)
#' var_of_interest <- survival_event[case_tag == "T"]
#' plot_mapper_with_pie_chart(mapper_object, var_of_interest, "Exitus")}
plot_mapper_with_pie_chart <- function(mapper_object, var_of_interest,
                                       name_var_of_interest = "Variable of interest",
                                       show_node_names = FALSE) {

  check_plot_vectors(mapper_object, var_of_interest, "character")
  # Create contingency table
  table_node_var <- percent_table_nodes(mapper_object, var_of_interest)

  adj_matrix <- mapper_object$adj_matrix
  # Remove diagonal from '1s'
  for (i in 1:dim(adj_matrix)[1]){
    adj_matrix[i,i] <- 0
  }

  g <- igraph::graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
  colors <- list(RColorBrewer::brewer.pal(length(unique(var_of_interest)),
                                          "Set3"))
  if (!show_node_names){node_names <- NA}
  else {node_names <- igraph::V(g)$name}

  plot(g, vertex.shape = "pie",
       vertex.pie = split(table_node_var, row(table_node_var)),
       vertex.pie.color = colors,
       vertex.size = 10,
       vertex.label = node_names,
       vertex.label.cex = 0.5,
       vertex.label.family = "sans",
       vertex.label.color = "black",
       vertex.label.dist = 2)

  graphics::legend(x = -2.5, y = 0.5,
                   legend = colnames(table_node_var),
                   fill = colors[[1]],
                   title = name_var_of_interest,
                   border = "black",
                   bty = "n",
                   title.cex = 0.7,
                   cex = 0.5)
}

