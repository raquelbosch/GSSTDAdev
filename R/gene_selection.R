#' @title Gene selection for partial clustering in Mapper based on variability.
#' @description
#' It selects genes for mapper based on different measures of variability.
#' @param case_disease_component Disease component matrix (output of the function
#' \code{generate_disease_component}) having selected only the columns
#' belonging to disease samples. The names of the rows must be the names of the
#' genes.
#' @param num_gen_select Number of genes to be selected (those with higher values
#' of the chosen variability measure).
#' @param gene_select_mapper_metric Measures of variability to choose as gene selection
#' criteria for Mapper. Choose as selection criteria between "mad" (mean
#' absolute deviation), "sd" (standard deviation), “iqr” (interquartile range)
#' and "mean_sd" (to choose genes with high mean and standard deviation
#' simultaneously).
#' @return Character vector with the names of the selected genes.
#' @import matrixStats
#' @import stats
gene_selection_mapper_by_variability <- function(case_disease_component,
                                                 num_gen_select,
                                                 gene_select_mapper_metric){
  if (gene_select_mapper_metric == "mad"){
    measure <- apply(case_disease_component, 1, stats::mad)
  }
  if (gene_select_mapper_metric == "sd"){
    measure <- apply(case_disease_component, 1, stats::sd)
  }
  if (gene_select_mapper_metric == "iqr"){
    measure <- matrixStats::rowIQRs(case_disease_component)
  }
  if (gene_select_mapper_metric == "mean_sd"){
    row_mean <- apply(case_disease_component, 1, mean)
    row_mean_scaled <- (row_mean - min(row_mean)) / (max(row_mean)
                                                     - min(row_mean))
    row_sd <- apply(case_disease_component, 1, stats::sd)
    row_sd_scaled <- (row_sd - min(row_sd)) / (max(row_sd) - min(row_sd))
    measure <- row_mean + row_sd
  }
  genes_selected <- names(measure[order(measure,decreasing = T)])[1:num_gen_select]

  return(genes_selected)
}


#' @title Gene selection for partial clustering in Mapper based on variability
#' and the relationship to survival.
#' @description
#' It selects genes for mapper based on the product of standard deviation
#' of the rows (genes) in the disease component matrix
#' plus one times the Z score obtained by fitting a cox proportional
#' hazard model to the level of each gene.
#' @param case_disease_component Disease component matrix (output of the function
#' \code{generate_disease_component}) having selected only the columns
#' belonging to disease samples. The names of the rows must be the names of the genes.
#' @param cox_all_matrix Output from the \code{cox_all_genes} function. Data.frame with
#' information on the relationship between genes and survival.
#' @param gene_select_surv_type Option. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis).
#' @param num_gen_select Number of genes to be selected (those with the highest
#' product value).
#' @return Character vector with the names of the selected genes.
gene_selection_by_sd_surv <- function(case_disease_component, cox_all_matrix,
                                      gene_select_surv_type, num_gen_select){
  # Same operation to both methods
  probes_test <- apply(case_disease_component, 1,stats::sd)+1

  cox_vector <- cox_all_matrix[,"z"]
  cox_vector[cox_vector < 0] <- ifelse(!is.na(cox_vector[cox_vector < 0]),
                                       cox_vector[cox_vector < 0] - 1,
                                       cox_vector[cox_vector < 0])
  cox_vector[cox_vector >= 0] <- ifelse(!is.na(cox_vector[cox_vector >= 0]),
                                        cox_vector[cox_vector >= 0] + 1,
                                        cox_vector[cox_vector >= 0])

  if(gene_select_surv_type == "abs"){
    probes_test <- probes_test * abs(cox_vector)
    genes_selected <- names(probes_test[order(probes_test,decreasing = T)])[1:num_gen_select]
  }else{
    probes_test <- probes_test * cox_vector
    if(num_gen_select %% 2 == 0){ num_gen_select <- num_gen_select/2}
    else{ num_gen_select <- (num_gen_select + 1)/2}

    genes_selected <- names(c(probes_test[order(probes_test,decreasing = T)][1:num_gen_select],
                              probes_test[order(probes_test,decreasing = F)][1:num_gen_select]))
  }
  return(genes_selected)
}




