#' @title Gene Structure Survival using Topological Data Analysis (GSSTDA).
#'
#' @description Gene Structure Survival using Topological Data Analysis.
#' This function implements an analysis for expression array data
#' based on the *Progression Analysis of Disease* developed by Nicolau
#' *et al.* (doi: 10.1073/pnas.1102826108) that allows the information
#' contained in an expression matrix to be condensed into a combinatory graph.
#' The novelty is that information on survival is integrated into the analysis.
#'
#' The analysis consists of 3 parts: a preprocessing of the data, the gene
#' selection and the filter function, and the mapper algorithm. The
#' preprocessing is specifically the Disease Specific Genomic Analysis (proposed
#' by Nicolau *et al.*) that consists of, through linear models, eliminating the
#' part of the data that is considered "healthy" and keeping only the component
#' that is due to the disease. Subsequently, the filter function is calculated,
#' whose values capture the survival associated with each patient, in addition
#' to selecting those genes that will be used in the Mapper algorithm.
#' Finally, the Mapper algorithm is applied from the disease component matrix of
#' pathological samples and the values of the filter function obtaining a
#' combinatory graph.
#' @param full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' @param survival_time Numerical vector of the same length as the number of
#' columns of \code{full_data}. In addition, the patients must be in the same
#' order as in \code{full_data}. For the patients whose sample is pathological
#' should be indicated the time between the disease diagnosis and event
#' (death, relapse or other). If the event has not occurred, it should be
#' indicated the time until the end of follow-up. Patients whose sample is
#' from healthy tissue must have an NA value
#' @param survival_event Numerical vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. For the the patients with pathological sample should
#' be indicated whether the event has occurred (1) or not (0). Only these
#' values are valid and healthy patients must have an NA value.
#' @param case_tag Character vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. It must be indicated for each patient whether its
#' sample is from pathological or healthy tissue. One value should be used to
#' indicate whether the patient's sample is healthy and another value should
#' be used to indicate whether the patient's sample is pathological.
#' The user will then be asked which one indicates whether the patient is
#' healthy. Only two values are valid in the vector in total.
#' @param control_tag Tag of the healthy sample.E.g. "T"
#' @param gamma A parameter that indicates the magnitude of the noise assumed in
#' the flat data matrix for the generation of the Healthy State Model. If it
#' takes the value `NA` the magnitude of the noise is assumed to be unknown.
#' By default gamma is unknown.
#' @param gene_select_surv_type Option. Options on how to select the genes to be
#' used in the calculation of the values of the filter function (and in the gene
#' selection for Mapper if the option "sd_surv" has been chosen) Select the
#' "Abs" option, which means that the genes with the highest absolute value are
#' chosen, or the "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis). "Top_Bot" default option.
#' @param percent_gen_select_for_fun_filt Percentage (from zero to one hundred)
#' of genes to be selected to be used in the calculation of the values of the
#' filter function. 1 default option.
#' @param gene_select_mapper_metric Gene selection criteria for Mapper. Choose as
#' selection criteria between:
#' - "mad": those genes in the disease component matrix with the highest mean
#' absolute deviation will be selected
#' - "sd": those genes with the highest standard deviation will be selected
#' - “iqr": those genes with the highest interquartile range will be selected
#' - "mean_sd": to choose genes with high mean and standard deviation simultaneously
#' - "sd_surv": these option selects genes for mapper based on the product of standard
#' deviation of the genes in the disease component matrix plus one times
#' the Z score obtained by fitting a cox proportional hazard model to the level
#' of each gene. "mad" default option.
#' @param percent_gen_select_for_mapper Percentage (from zero to one hundred)
#' of genes to be selected to be used in Mapper. 5 default option.
#' @param num_intervals Parameter for the mapper algorithm. Number of
#' intervals used to create the first sample partition based on
#' filtering values. By default the root of the number of individuals included
#' as input to Mapper
#' @param percent_overlap Parameter for the mapper algorithm. Percentage
#' of overlap between intervals. Expressed as a percentage. 40 default option.
#' @param type_covering Parameter to choose how to construct the covering.
#' Choose between "classical" or "uniform". If the ‘classic’ option is
#' selected, the covering is formed by splitting the range of filter function
#' values into overlapping intervals of the same length. If the ‘uniform’ option
#' is selected, the covering will be formed with overlapping intervals of the
#' filter function values containing approximately the same number of
#' individuals. "uniform" default option.
#' @param distance_type Parameter for the mapper algorithm.
#' Type of distance to be used for clustering. Choose between correlation
#' ("correlation") and euclidean ("euclidean"). "correlation" default option.
#' @param clustering_type Parameter for the mapper algorithm. Type of
#' clustering method. Choose between "hierarchical" and "PAM"
#' (“partition around medoids”) options. "hierarchical" default option.
#' @param num_bins_when_clustering Parameter for the mapper algorithm.
#' Number of bins to generate the histogram employed by the standard
#' optimal number of cluster finder method. Parameter not necessary if the
#' "optimal_clust_mode" option is "silhouette" or the "clust_type" is "PAM".
#' 8 default option.
#' @param linkage_type Parameter for the mapper algorithm. Linkage criteria
#' used in hierarchical clustering. Choose between "single" for single-linkage
#' clustering, "complete" for complete-linkage clustering, "average" for
#' average linkage clustering (or UPGMA) or "ward.D" for Ward's method. Only
#' necessary for hierarchical clustering. "ward.D" default option.
#' @param optimal_clustering_mode Method for selection optimal number of
#' clusters. It is only necessary if the chosen type of algorithm is
#' hierarchical. In this case, choose between "standard" (the method used
#' in the original mapper article) or "silhouette". In the case of the PAM
#' algorithm, the method will always be "silhouette".
#' @param silhouette_threshold Minimum value of \eqn{\overline{s}}{s-bar} that a set of
#' clusters must have to be chosen as optimal. Within each interval of the
#' filter function, the average silhouette values \eqn{\overline{s}}{s-bar} are computed
#' for all possible partitions from $2$ to $n-1$, where $n$ is the number of
#' samples within a specific interval. The $n$ that produces the highest value
#' of \eqn{\overline{s}}{s-bar} and that exceeds a specific threshold is selected as the
#' optimum number of clusters. If no partition produces an \eqn{\overline{s}}{s-bar}
#' exceeding the chosen threshold, all samples are then assigned to a unique
#' cluster. The default value is $0.25$. The threshold of $0.25$ for
#' \eqn{\overline{s}}{s-bar} has been chosen based on standard practice, recognizing it
#' as a moderate value that reflects adequate separation and cohesion within
#' clusters.
#' @param dim_reduction Boolean. Parameter indicating whether or not centering,
#' scaling and dimensionality reduction of the data is wanted prior to partial
#' clustering in Mapper. In dimensionality reduction, principal components that
#' explain at least 80% of the variance are chosen. This process is performed on
#' each subset composed by the columns of the individuals that are part of each
#' interval. TRUE indicates that centering, scaling and dimensionality reduction
#' will be performed. FALSE means that this step will be omitted and clustering
#' will be performed without transforming the data.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{gsstda} object. It contains:
#' - the matrix with the normal space \code{normal_space},
#' - the matrix of the disease components normal_space \code{matrix_disease_component},
#' - a matrix with the results of the application of proportional hazard models
#' for each gene (\code{cox_all_matrix)},
#' - the selected genes for mapper (\code{genes_selected_mapper}) and for the
#' calculation of the values of the filter function (\code{genes_selected_fun_filt}),
#' - \code{case_genes_disease_component} (the matrix of disease components with
#' only the rows of the selected genes and only the columns of the pathological
#' samples which is one of the inputs for the \code{mapper} function),
#' - and a \code{mapper_obj} object. This \code{mapper_obj} object contains the
#' values of the intervals (interval_data), the samples included in each
#' interval (sample_in_level), information about the cluster to which the
#' individuals in each interval belong (clustering_all_levels), a list including
#' the individuals contained in each detected node (node_samples), a list
#' including the nodes contained in each interval (nodes_in_lev), their size
#' (node_sizes), the average of the filter function values of the individuals
#' of each node (node_average_filt) and the adjacency matrix linking the nodes
#' (adj_matrix). Moreover, information is provided on the number of nodes,
#' the average node size, the standard deviation of the node size, the number
#' of connections between nodes, the proportion of connections to all possible
#' connections and the number of ramifications.
#' @export
#' @examples
#' \donttest{
#' gsstda_object <- gsstda(full_data,  survival_time, survival_event, case_tag, gamma=NA,
#'                  gene_select_surv_type="Top_Bot", percent_gen_select_for_fun_filt=1,
#'                  gene_select_mapper_metric="mad", percent_gen_select_for_mapper=5,
#'                  num_intervals = 4, percent_overlap = 50,
#'                  type_covering = "uniform", distance_type = "euclidean",
#'                  num_bins_when_clustering = 8, clustering_type = "hierarchical",
#'                  linkage_type = "single")}
gsstda <- function(full_data, survival_time, survival_event, case_tag,
                   control_tag=NA, gamma=NA, gene_select_surv_type="Top_Bot",
                   percent_gen_select_for_fun_filt=1, gene_select_mapper_metric="mad",
                   percent_gen_select_for_mapper=5, num_intervals=NULL, percent_overlap=40,
                   type_covering = "uniform", distance_type="correlation",
                   clustering_type="hierarchical", num_bins_when_clustering = 8,
                   linkage_type = "ward.D", optimal_clustering_mode = NA,
                   silhouette_threshold = 0.25, dim_reduction=TRUE, na.rm=TRUE){
  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  full_data <- check_full_data(full_data, na.rm)
  #Select the control_tag. This do it inside of the dsga function
  #Check and obtain gene selection (we use in the gene_select_surv). It execute in Block II
  #return_check_gene_selection <- check_gene_selection(nrow(full_data),
  #                                                    gene_select_surv_type,
  #                                                    percent_gen_select_for_fun_filt,
  #                                                    gene_select_mapper_metric,
  #                                                    percent_gen_select_for_mapper)
  #num_gen_select_for_fun_filt <- return_check_gene_selection[[1]]
  #num_gen_select_for_mapper <- return_check_gene_selection[[2]]

  #Don't check filter_values because it is not created.
  filter_values <- c()

  check_return <- check_arg_mapper(full_data, filter_values, type_covering,
                                   distance_type, clustering_type,
                                   linkage_type, optimal_clustering_mode,
                                   silhouette_threshold, dim_reduction, na.rm)

  full_data <- check_return[[1]]
  filter_values <- check_return[[2]]
  optimal_clustering_mode <- check_return[[3]]


  ################### BLOCK I: Pre-process. dsga (using "NT" control_tag) ##############################
  dsga_obj <- dsga(full_data, survival_time, survival_event, case_tag, control_tag, gamma, na.rm = "checked")
  matrix_disease_component <- dsga_obj[["matrix_disease_component"]]
  control_tag <- dsga_obj[["control_tag"]]
  full_data <- dsga_obj[["full_data"]]
  survival_event <- dsga_obj[["survival_event"]]
  survival_time <- dsga_obj[["survival_time"]]
  case_tag <- dsga_obj[["case_tag"]]

  ################### BLOCK II: Gene selection (using "T" control_tag) ##################################
  gene_selection_object <- gene_selection(dsga_obj, gene_select_surv_type,
                                          percent_gen_select_for_fun_filt,
                                          gene_select_mapper_metric,
                                          percent_gen_select_for_mapper)
  cox_all_matrix <- gene_selection_object[["cox_all_matrix"]]
  genes_selected_for_mapper <- gene_selection_object[["genes_selected_for_mapper"]]
  genes_selected_for_fun_filt <- gene_selection_object[["genes_selected_for_fun_filt"]]
  case_genes_disease_component <- gene_selection_object[["case_genes_disease_component"]]
  filter_values <- gene_selection_object[["filter_values"]]

  ################### BLOCK III: Create mapper object where the arguments are checked ###################
  message("\nBLOCK III: The mapper process is started")

  # Transpose case_genes_disease_component: rows = patient, columns = genes
  #case_genes_disease_component <- t(case_genes_disease_component)

  #   Check filter_values
  check_filter <- check_filter_values(case_genes_disease_component, filter_values)
  case_genes_disease_component <- check_filter[[1]]
  filter_values <- check_filter[[2]]

  mapper_obj <- mapper(case_genes_disease_component, filter_values,
                       num_intervals, percent_overlap, type_covering,
                       distance_type, clustering_type, num_bins_when_clustering,
                       linkage_type, optimal_clustering_mode,
                       silhouette_threshold, dim_reduction = dim_reduction,
                       na.rm = "checked")

  message("\nBLOCK III: The mapper process is finished")


  ############################################  Create the object #########################################
  gsstda_object <- list("normal_space" = dsga_obj[["normal_space"]],
                        "matrix_disease_component" = matrix_disease_component,
                        "cox_all_matrix" = cox_all_matrix,
                        "genes_selected_for_mapper" = genes_selected_for_mapper,
                        "genes_selected_for_fun_filt" = genes_selected_for_fun_filt,
                        "case_genes_disease_component" = case_genes_disease_component,
                        "mapper_obj" = mapper_obj
                        )

  class(gsstda_object) <- "gsstda_obj"
  return(gsstda_object)
}
