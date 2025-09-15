#' @title Disease-Specific Genomic Analysis
#' @description Disease-Specific Genomic Analysis (dsga).
#' This analysis, developed by Nicolau *et al.*, allows the calculation of
#' the "disease component" of a expression matrix which consists of, through
#' linear models, eliminating the part of the data  that is considered normal
#' or healthy and keeping only the component that is due to the disease. It
#' is intended to precede other techniques like classification or clustering.
#' For more information see *Disease-specific genomic analysis: identifying
#' the signature of pathologic biology* (doi: 10.1093/bioinformatics/btm033).
#' @param full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' @param case_tag Character vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. It must be indicated for each patient whether its
#' sample is from pathological or healthy tissue. One value should be used to
#' indicate whether the patient's sample is healthy and another value should
#' be used to indicate whether the patient's sample is pathological.
#' The user will then be asked which one indicates whether the patient is
#' healthy. Only two values are valid in the vector in total.
#' @param control_tag Tag of the healthy sample. E.g. "T"
#' @param gamma A parameter that indicates the magnitude of the noise assumed in
#' the flat data matrix for the generation of the Healthy State Model. If it
#' takes the value `NA` the magnitude of the noise is assumed to be unknown.
#' By default gamma is unknown.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{dsga} object. It contains: the \code{full_data} without
#' NAN's values, the label designated for healthy samples (\code{control_tag}),
#' the \code{case_tag} vector without NAN's values, the matrix with the normal
#' space (linear space generated from normal tissue samples) and the matrix of
#' the disease components (the transformed full_data matrix from which the
#' normal component has been removed).
#' @export
#' @examples
#' \donttest{
#' dsga_obj <- dsga(full_data, case_tag)}
dsga <- function(full_data, case_tag, control_tag = NA, gamma = NA,
                 na.rm = TRUE){
  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  if (na.rm != "checked"){
    full_data <- check_full_data(full_data, na.rm)

    #Select the control_tag
    return_check <- check_case_tag(full_data, case_tag, control_tag)
    control_tag <- return_check[[1]]
    full_data <- return_check[[2]]
    case_tag <- return_check[[3]]
    }

  ################### BLOCK I: Pre-process. dsga (using "NT" control_tag) ##############################
  message("\nBLOCK I: The pre-process dsga is started")
  #   Select the normal tissue data gene expression matrix.
  normal_tiss <- full_data[,which(case_tag == control_tag)]

  #   Obtain the gene expression matrix containing the flattened version of the vectors.
  matrix_flatten_normal_tiss <- flatten_normal_tiss(normal_tiss)
  #   Obtain the normal space
  normal_space <- denoise_rectangular_matrix(matrix_flatten_normal_tiss, gamma)
  #   Obtain the disease component of the normal_space
  matrix_disease_component <- generate_disease_component(full_data, normal_space)

  message("\nBLOCK I: The pre-process dsga is finished\n")

  ############################################  Create the object #########################################
  dsga_object <- list("full_data" = full_data,
                      "control_tag" = control_tag,
                      "case_tag" = case_tag,
                      "normal_space" = normal_space,
                      "matrix_disease_component" = matrix_disease_component)

  class(dsga_object) <- "dsga_object"

  return(dsga_object)
}


#' @title Gene selection and filter function
#' @description Gene selection and calculation of filter function values.
#' After fitting a Cox proportional hazard model to each gene, this function
#' makes a selection of genes according to their relationship with survival.
#' Subsequently, with the genes selected, the values of the filtering functions
#' are calculated for each patient. The filter function allows to summarise each
#' vector of each individual in a single data. Specifically, it computes a
#' multivariate Cox proportional hazard model using all previously selected
#' genes as predictor variables. Subsequently, it calculates for each patient
#' using this model their linear predictor, which corresponds to the exponent
#' of their proportionality constant.
#' Finally, this function selects the genes,
#' i.e. the columns, from the disease component matrix that will be used as
#' input data for the mapper algorithm based on the selected criteria.
#' @param data_object A dsga_object or an object with:
#' - full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' - case_tag Character vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. It must be indicated for each patient whether its
#' sample is from pathological or healthy tissue. One value should be used to
#' indicate whether the patient's sample is healthy and another value should
#' be used to indicate whether the patient's sample is pathological.
#' The user will then be asked which one indicates whether the patient is
#' healthy. Only two values are valid in the vector in total.
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
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{gene_selection_object}. It contains:
#' - the matrix with which the gene selection has been performed without NAN's
#' values (\code{data}). It is the \code{matrix_disease_component} in case it has been
#' performed from a \code{dsga_object} or \code{full_data} in the opposite case.
#' - the \code{cox_all_matrix} (a matrix with the results of the application of
#' proportional hazard models: with the regression coefficients, the odds ratios,
#' the standard errors of each coefficient, the Z values (coef/se_coef) and
#' the p-values for each Z value)
#' - a vector with the name of the selected genes for mapper
#' (\code{genes_selected_mapper}) and for the calculation of the values of the
#' filter function (\code{genes_selected_fun_filt}).
#' - \code{case_genes_disease_component}: the matrix of disease components with
#' only the rows of the selected genes and only the columns of the pathological
#' (case) samples. This matrix is one of the inputs for the \code{mapper} function.
#' - and the vector of the values of the filter function of pathological samples
#' \code{filter_values}.
#' @export
#' @examples
#' \donttest{
#' data_object <- list("full_data" = full_data, "case_tag" = case_tag)
#' class(data_object) <- "data_object"
#' gene_selection_obj <- gene_selection(data_object, survival_time,
#' survival_event, gene_select_surv_type ="top_bot",
#' percent_gen_select_for_fun_filt=1, gene_select_mapper_metric="mad",
#' percent_gen_select_for_mapper=5)}
gene_selection <- function(data_object, survival_time, survival_event,
                           gene_select_surv_type = "Top_Bot",
                           percent_gen_select_for_fun_filt = 1,
                           gene_select_mapper_metric = "mad",
                           percent_gen_select_for_mapper = 5,
                           na.rm = TRUE){
  UseMethod("gene_selection")
}

#' @title Private gene_selection_
#' @description Private function to gene selection
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
#' @param control_tag_cases Numeric vector with the indices of the columns of
#' \code{full_data} and/or \code{matrix_disease_component}
#' corresponding to the healthy sample patients.
#' @param gene_select_surv_type Option. Options on how to select the genes to be
#' used in the calculation of the values of the filter function (and in the gene
#' selection for Mapper if the option "sd_surv" has been chosen). Select the
#' "Abs" option, which means that the genes with the highest absolute value are
#' chosen, or the "Top_Bot" option, which means that half of the selected genes
#' are those with the highest value (positive value, i.e. worst survival
#' prognosis) and the other half are those with the lowest value (negative
#' value, i.e. best prognosis). "Top_Bot" default option.
#' @param num_gen_select_for_fun_filt Number of genes to be selected to be used
#' for the calculation of the values of the filter function.
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
#' of each gene.
#' @param num_gen_select_for_mapper Number of genes to be selected to be used
#' in mapper.
#' @param matrix_disease_component Optional, only necessary in case of gene
#' selection after DSGA has been performed. Matrix of the disease components
#' (the transformed \code{full_data} matrix from which the normal component has
#' been removed) from the \code{dsga_function}.
#' @return A \code{gene_selection_object}. It contains:
#' - the matrix with which the gene selection has been performed without NAN's
#' values (\code{data}). It is the \code{matrix_disease_component} in case it has been
#' performed from a \code{dsga_object} or \code{full_data} in the opposite case.
#' - the \code{cox_all_matrix} (a matrix with the results of the application of
#' proportional hazard models: with the regression coefficients, the odds ratios,
#' the standard errors of each coefficient, the Z values (coef/se_coef) and
#' the p-values for each Z value)
#' - a vector with the name of the selected genes for mapper
#' (\code{genes_selected_mapper}) and for the calculation of the values of the
#' filter function (\code{genes_selected_fun_filt}).
#' - \code{case_genes_disease_component}: the matrix of disease components with
#' only the rows of the selected genes and only the columns of the pathological
#' (case) samples. This matrix is one of the inputs for the \code{mapper} function.
#' - and the vector of the values of the filter function of pathological samples
#' \code{filter_values}.
#' @export
#' @examples
#' \donttest{
#' control_tag_cases <- which(case_tag == "NT")
#' gene_selection_obj <- gene_selection_(full_data, survival_time, survival_event,
#' control_tag_cases, gene_select_surv_type ="top_bot", num_gen_select_for_fun_filt = 200,
#' gene_select_mapper_metric="mad", num_gen_select_for_mapper = 1000)}
gene_selection_ <- function(full_data, survival_time, survival_event,
                            control_tag_cases, gene_select_surv_type,
                            num_gen_select_for_fun_filt,
                            gene_select_mapper_metric,
                            num_gen_select_for_mapper,
                            matrix_disease_component = NULL){

  message("\nBLOCK II: The gene selection is started\n")

  if(is.null(matrix_disease_component)) {
    matrix_disease_component <- full_data
  }

  #Remove NAN's values (case_tag == control_tag) of survival_time and survival_event
  survival_time <- survival_time[-control_tag_cases]
  survival_event <- survival_event[-control_tag_cases]
  #Select the disease component of the "T" control_tag
  case_full_data <- full_data[,-control_tag_cases]
  case_disease_component <- matrix_disease_component[,-control_tag_cases]

  # Univariate cox proportional hazard models for the expression levels of each gene included in the
  #provided dataset
  cox_all_matrix <- cox_all_genes(case_full_data, survival_time, survival_event)

  # Select genes in full_data for filter function calulation
  genes_selected_fun_filt <- gene_selection_fun_filt(cox_all_matrix,
                                                     gene_select_surv_type,
                                                     num_gen_select_for_fun_filt)
  # Selects genes for mapper
  if (gene_select_mapper_metric == "sd_surv"){
    genes_selected_mapper <- gene_selection_by_sd_surv(case_disease_component,
                                                       cox_all_matrix, gene_select_surv_type,
                                                       num_gen_select_for_mapper)
  } else {
    genes_selected_mapper <- gene_selection_mapper_by_variability(case_disease_component,
                                                                  num_gen_select_for_mapper,
                                                                  gene_select_mapper_metric)
  }

  # Calculate the values of the filter function, it is necessary to select
  # the rows of the genes chosen
  filter_values <- calculate_filter_function(case_full_data[genes_selected_fun_filt,],
                                             survival_time, survival_event)

  message("\nBLOCK II: The gene selection is finished\n")

  # Select genes in matrix_disease_component or full_data (if don't apply Block I)
  # for mapper
  case_genes_disease_component <- matrix_disease_component[genes_selected_mapper,
                                                           -control_tag_cases]

  gene_selection_object <- list("data" = matrix_disease_component,
                                "cox_all_matrix" = cox_all_matrix,
                                "genes_selected_for_mapper" = genes_selected_mapper,
                                "genes_selected_for_fun_filt" = genes_selected_fun_filt,
                                "case_genes_disease_component" = case_genes_disease_component,
                                "filter_values" = filter_values
  )
  class(gene_selection_object) <- "gene_selection_object"

  return(gene_selection_object)
}

#' @title gene_selection_classes.dsga_object
#' @description Private function to select Gene with dsga object
#' @param data_object dsga object information
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
#' @param gene_select_surv_type Option. Options on how to select the genes to be
#' used in the mapper. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
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
#' - "iqr": those genes with the highest interquartile range will be selected
#' - "mean_sd": to choose genes with high mean and standard deviation simultaneously
#' - "sd_surv": these option selects genes for mapper based on the product of standard
#' deviation of the genes in the disease component matrix plus one times
#' the Z score obtained by fitting a cox proportional hazard model to the level
#' of each gene.
#' @param percent_gen_select_for_mapper Percentage (from zero to one hundred)
#' of genes to be selected to be used in Mapper. 5 default option.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{gene_selection_object}. It contains:
#' - the matrix with which the gene selection has been performed without NAN's
#' values (\code{data}). It is the \code{matrix_disease_component} in case it has been
#' performed from a \code{dsga_object} or \code{full_data} in the opposite case.
#' - the \code{cox_all_matrix} (a matrix with the results of the application of
#' proportional hazard models: with the regression coefficients, the odds ratios,
#' the standard errors of each coefficient, the Z values (coef/se_coef) and
#' the p-values for each Z value)
#' - a vector with the name of the selected genes for mapper
#' (\code{genes_selected_mapper}) and for the calculation of the values of the
#' filter function (\code{genes_selected_fun_filt}).
#' - \code{case_genes_disease_component}: the matrix of disease components with
#' only the rows of the selected genes and only the columns of the pathological
#' (case) samples. This matrix is one of the inputs for the \code{mapper} function.
#' - and the vector of the values of the filter function of pathological samples
#' \code{filter_values}.
#' @export
#' @examples
#' \donttest{
#' dsga_obj <- dsga(full_data, case_tag)
#' gene_selection_object <- gene_selection(dsga_obj, survival_time,
#'                                         survival_event,
#'                                         gene_select_surv_type ="top_bot",
#'                                         percent_gen_select_for_fun_filt=1,
#'                                         gene_select_mapper_metric="mad",
#'                                         percent_gen_select_for_mapper=5)}
gene_selection.dsga_object <- function(data_object, survival_time,
                                       survival_event, gene_select_surv_type,
                                       percent_gen_select_for_fun_filt,
                                       gene_select_mapper_metric,
                                       percent_gen_select_for_mapper,
                                       na.rm = TRUE){
  print(class(data_object))

  matrix_disease_component <- data_object[["matrix_disease_component"]]
  full_data <- data_object[["full_data"]]

  if (na.rm != "checked"){
    check_surv_vectors(full_data, survival_time, survival_event)
  }

  # Check and obtain gene selection. This function is not executed within the
  # gsstda function
  return_check_gene_selection <- check_gene_selection(nrow(matrix_disease_component),
                                                      gene_select_surv_type,
                                                      percent_gen_select_for_fun_filt,
                                                      gene_select_mapper_metric,
                                                      percent_gen_select_for_mapper)
  num_gen_select_for_fun_filt <- return_check_gene_selection[[1]]
  num_gen_select_for_mapper <- return_check_gene_selection[[2]]

  control_tag <- data_object[["control_tag"]]
  case_tag <- data_object[["case_tag"]]

  control_tag_cases <- which(case_tag == control_tag)
  gene_selection_object <- gene_selection_(full_data, survival_time,
                                           survival_event, control_tag_cases,
                                           gene_select_surv_type, num_gen_select_for_fun_filt,
                                           gene_select_mapper_metric, num_gen_select_for_mapper,
                                           matrix_disease_component)

  return(gene_selection_object)
}

#' @title gene_selection_classes.default
#' @description Private function to select Gene without dsga process
#' @param data_object Object with:
#' - full_data Input matrix whose columns correspond to the patients and
#' rows to the genes.
#' - case_tag Character vector of the same length as the number of
#' columns of \code{full_data}. Patients must be in the same order as in
#' \code{full_data}. It must be indicated for each patient whether its
#' sample is from pathological or healthy tissue. One value should be used to
#' indicate whether the patient's sample is healthy and another value should
#' be used to indicate whether the patient's sample is pathological.
#' The user will then be asked which one indicates whether the patient is
#' healthy. Only two values are valid in the vector in total.
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
#' @param gene_select_surv_type Option. Options on how to select the genes to be
#' used in the mapper. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
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
#' of each gene.
#' @param percent_gen_select_for_mapper Percentage (from zero to one hundred)
#' of genes to be selected to be used in Mapper. 5 default option.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{gene_selection_object}. It contains:
#' - the matrix with which the gene selection has been performed without NAN's
#' values (\code{data}). It is the \code{matrix_disease_component} in case it has been
#' performed from a \code{dsga_object} or \code{full_data} in the opposite case.
#' - the \code{cox_all_matrix} (a matrix with the results of the application of
#' proportional hazard models: with the regression coefficients, the odds ratios,
#' the standard errors of each coefficient, the Z values (coef/se_coef) and
#' the p-values for each Z value)
#' - a vector with the name of the selected genes for mapper
#' (\code{genes_selected_mapper}) and for the calculation of the values of the
#' filter function (\code{genes_selected_fun_filt}).
#' - \code{case_genes_disease_component}: the matrix of disease components with
#' only the rows of the selected genes and only the columns of the pathological
#' (case) samples. This matrix is one of the inputs for the \code{mapper} function.
#' - and the vector of the values of the filter function of pathological samples
#' \code{filter_values}.
#' @export
#' @examples
#' \donttest{
#' data_object <- list("full_data" = full_data, "case_tag" = case_tag)
#' class(data_object) <- "data_object"
#' gene_selection_obj <- gene_selection(data_object, survival_time,
#'                                      survival_event,
#'                                      gene_select_surv_type ="top_bot",
#'                                      percent_gen_select_for_fun_filt=1,
#'                                      gene_select_mapper_metric="mad",
#'                                      percent_gen_select_for_mapper=5)}
gene_selection.default <- function(data_object, survival_time,
                                   survival_event, gene_select_surv_type,
                                   percent_gen_select_for_fun_filt,
                                   gene_select_mapper_metric,
                                   percent_gen_select_for_mapper, na.rm = TRUE){
  full_data <- data_object[["full_data"]]
  case_tag <- data_object[["case_tag"]]

  ################################ Prepare data and check data ########################################
  #Check the arguments introduces in the function
  full_data <- check_full_data(full_data, na.rm)

  #Check and obtain gene selection
  return_check_gene_selection <- check_gene_selection(nrow(full_data),
                                                      gene_select_surv_type,
                                                      percent_gen_select_for_fun_filt,
                                                      gene_select_mapper_metric,
                                                      percent_gen_select_for_mapper)
  num_gen_select_for_fun_filt <- return_check_gene_selection[[1]]
  num_gen_select_for_mapper <- return_check_gene_selection[[2]]

  #Select the control_tag
  return_check <- check_vectors(full_data, survival_time, survival_event,
                                case_tag, control_tag = NA)

  control_tag <- return_check[[1]]
  full_data <- return_check[[2]]
  survival_time <- return_check[[3]]
  survival_event <- return_check[[4]]
  case_tag <- return_check[[5]]

  control_tag_cases <- which(case_tag == control_tag)

  gene_selection_object <- gene_selection_(full_data, survival_time,
                                           survival_event, control_tag_cases,
                                           gene_select_surv_type, num_gen_select_for_fun_filt,
                                           gene_select_mapper_metric, num_gen_select_for_mapper)

  return(gene_selection_object)
}


#' @title Mapper object
#'
#' @description TDA are persistent homology and mapper. Persistent homology
#' borrows ideas from abstract algebra to identify particular aspects
#' related to the shape of the data such as the number of connected
#' components and the presence of higher-dimensional holes, whereas
#' mapper condenses the information of high-dimensional datasets into
#' a combinatory graph or simplicial complex that is referred to as
#' the skeleton of the dataset. This implementation is the mapper of one
#' dimension, i.e. using only one filter function value.
#' @param data Input matrix whose columns correspond to the individuals
#' and rows to the features.
#' @param filter_values Vector obtained after applying the filtering function
#' to the input matrix, i.e, a vector with the filtering function
#' values for each included sample.
#' @param num_intervals Number of intervals used to create the first sample
#' partition based on filtering values. By default the root of the number of
#' individuals included as input to Mapper.
#' @param percent_overlap Percentage of overlap between intervals. Expressed
#' as a percentage. 40 default option.
#' @param type_covering Parameter to choose how to construct the covering.
#' Choose between "classical" or "uniform". If the ‘classic’ option is
#' selected, the covering is formed by splitting the range of filter function
#' values into overlapping intervals of the same length. If the ‘uniform’ option
#' is selected, the covering will be formed with overlapping intervals of the
#' filter function values containing approximately the same number of
#' individuals. "uniform" default option.
#' @param distance_type Type of distance to be used for clustering.
#' Choose between correlation ("correlation") and euclidean ("euclidean"). "correlation"
#' default option.
#' @param clustering_type Type of clustering method. Choose between
#' "hierarchical" and "PAM" (“partition around medoids”) options.
#' "hierarchical" default option.
#' @param num_bins_when_clustering Number of bins to generate the
#' histogram employed by the standard optimal number of cluster finder
#' method. Parameter not necessary if the "optimal_clustering_mode" option
#' is "silhouette" or the "clustering_type" is "PAM". 8 default option.
#' @param linkage_type Linkage criteria used in hierarchical clustering.
#' Choose between "single" for single-linkage clustering, "complete" for
#' complete-linkage clustering, "average" for average linkage clustering
#' (or UPGMA) or "ward.D" for Ward's method. Only necessary for hierarchical
#' clustering. "ward.D" default option.
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
#' explain at least 80% of the variance are chosen. This process is performed
#' on each subset composed by the columns of the individuals that are part of
#' each interval. TRUE indicates that centering, scaling and dimensionality
#' reduction will be performed. FALSE means that this step will be omitted and
#' clustering will be performed without transforming the data.
#' @param na.rm \code{logical}. If \code{TRUE}, \code{NA} rows are omitted.
#' If \code{FALSE}, an error occurs in case of \code{NA} rows. TRUE default
#' option.
#' @return A \code{mapper_obj} object. It contains the values of the intervals
#' (interval_data), the samples included in each interval (sample_in_level),
#' information about the cluster to which the individuals in each interval
#' belong (clustering_all_levels), a list including the individuals contained
#' in each detected node (node_samples), a list including the nodes contained
#' in each interval (nodes_in_lev), their size (node_sizes), the
#' average of the filter function values of the individuals of each node
#' (node_average_filt) and the adjacency matrix linking the nodes (adj_matrix).
#' @export
#' @examples
#' \donttest{
#' control_tag_cases <- which(case_tag == "NT")
#' gene_selection_object <- gene_selection_(full_data, survival_time, survival_event,
#' control_tag_cases, gene_select_surv_type ="top_bot", num_gen_select_for_fun_filt = 200,
#' gene_select_mapper_metric="mad", num_gen_select_for_mapper = 1000)
#'
#' mapper_object <- mapper(data = gene_selection_object[["case_genes_disease_component"]],
#' filter_values = gene_selection_object[["filter_values"]],
#' num_intervals = 5,
#' percent_overlap = 40, type_covering = "uniform",
#' distance_type = "correlation", clustering_type = "hierarchical",
#' linkage_type = "single")}
mapper <- function(data, filter_values, num_intervals = NULL, percent_overlap = 40,
                   type_covering = "uniform", distance_type = "correlation",
                   clustering_type = "hierarchical", num_bins_when_clustering = 8,
                   linkage_type = "ward.D", optimal_clustering_mode = NA,
                   silhouette_threshold = 0.25, dim_reduction = TRUE,
                   na.rm = TRUE){
  # Don't call by GSSTDA function
  if (na.rm != "checked"){
    # Check the data introduces
    data <- check_full_data(data)
    # Check mapper arguments
    check_return <- check_arg_mapper(data, filter_values, type_covering,
                                     distance_type, clustering_type, linkage_type,
                                     optimal_clustering_mode, silhouette_threshold,
                                     dim_reduction)
    data <- check_return[[1]]
    filter_values <- check_return[[2]]
    optimal_clustering_mode <- check_return[[3]]
  }

  if(is.null(num_intervals)) {
    num_intervals <- round(sqrt(ncol(data)), digits = 0)
  }

  mapper_object_ini <- list("data" = data,
                            "filter_values" = filter_values,
                            "num_intervals" = num_intervals,
                            "percent_overlap" = percent_overlap/100,
                            "type_covering" = type_covering,
                            "distance_type" = distance_type,
                            "optimal_clustering_mode" = optimal_clustering_mode,
                            "num_bins_when_clustering" = num_bins_when_clustering,
                            "clustering_type" = clustering_type,
                            "linkage_type" = linkage_type,
                            "optimal_clustering_mode" = optimal_clustering_mode,
                            "silhouette_threshold" = silhouette_threshold,
                            "dim_reduction" = dim_reduction)

  class(mapper_object_ini) <- "mapper_initialization"

  mapper_object <- one_D_Mapper(mapper_object_ini)

  return(mapper_object)
}


