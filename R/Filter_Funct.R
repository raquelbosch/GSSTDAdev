#' @title Survival analysis based on gene expression levels.
#' @description It carries out univariate cox proportional hazard models for
#' the expression levels of each gene included in the provided dataset
#' (case_full_data) and their link with relapse-free or overall survival.
#' @param case_full_data Input matrix whose columns correspond to the patients
#' and rows to the genes, having selected only the columns belonging to disease
#' samples. The names of the rows must be the names of the genes.
#' @param survival_time Numeric vector that includes time to the event
#' information after removal of NAs corresponding to healthy tissue samples
#' (case_tag == control_tag).
#' @param survival_event Numeric vector that indicates if relapse or death
#' have been produced (0 and 1s) after removal of NAs corresponding to healthy
#' tissue samples (case_tag == control_tag)
#' @return A matrix with the results of the application of proportional
#' hazard models using the expression levels of each gene as covariate.
#' The \code{coef} column corresponds to the regression coefficient; the
#' \code{exp_coef} column corresponds to the value of e^coef  (which is
#' interpreted as the odds ratio); the \code{se_coef} column corresponds
#' to the standard error of each coefficient; the \code{Z} column corresponds
#' to the value of coef/se_coef (the higher the Z value, the higher the
#' significance of the variable) and the \code{Pr_z} column corresponds to
#' the p-value for each Z value.
#' @import survival
cox_all_genes <- function(case_full_data, survival_time, survival_event){
  message("Calculating the matrix of Zcox")
  pb <- utils::txtProgressBar(min = 0, max = nrow(case_full_data), style = 3)

  list_out <- list()
  for(i in 1:nrow(case_full_data)){
    utils::setTxtProgressBar(pb, i)

    temp <- summary(survival::coxph(survival::Surv(survival_time,survival_event)
                                    ~case_full_data[i,]))$coefficients[1,]
    list_out[[i]] <- temp
  }
  cox_all_matrix <- data.frame(do.call("rbind",list_out))
  colnames(cox_all_matrix) <-  c("coef","exp_coef","se_coef","z","Pr_z")
  rownames(cox_all_matrix) <- rownames(case_full_data)
  cox_all_matrix <- as.matrix(cox_all_matrix)
  return(cox_all_matrix)
}


#' @title Gene selection for filter function based on relationship to survival.
#' @description
#' It selects genes for mapper filter function based on Z score obtained by
#' fitting a cox proportional hazard model to the level of each gene.
#' @param cox_all_matrix Output from the \code{cox_all_genes} function. Data.frame with
#' information on the relationship between genes and survival.
#' @param gene_select_surv_type Option. Select the "Abs" option, which means that the
#' genes with the highest absolute value are chosen, or the
#' "Top_Bot" option, which means that half of the selected
#' genes are those with the highest value (positive value, i.e.
#' worst survival prognosis) and the other half are those with the
#' lowest value (negative value, i.e. best prognosis).
#' @param num_gen_select Number of genes to be selected (those with the highest
#' Z score).
#' @return Character vector with the names of the selected genes.
gene_selection_fun_filt <- function(cox_all_matrix,
                                    gene_select_surv_type, num_gen_select){
  z_score_vector <- cox_all_matrix[,"z"]

  if(gene_select_surv_type == "abs"){
    z_score_vector <- abs(z_score_vector)
    genes_selected <- names(z_score_vector[order(z_score_vector,decreasing = T)])[1:num_gen_select]
  }else{
    if(num_gen_select %% 2 == 0){ num_gen_select <- num_gen_select/2}
    else{ num_gen_select <- (num_gen_select + 1)/2}

    genes_selected <- names(c(z_score_vector[order(z_score_vector,decreasing = T)][1:num_gen_select],
                              z_score_vector[order(z_score_vector,decreasing = F)][1:num_gen_select]))
  }
  return(genes_selected)
}

#' @title Filtering function
#' @description
#' A filtering function for Mapper that projects $$R$^n$ into $R$.
#' It calculates for each column of the matrix (each patient), its value
#' of the filtering function. Specifically, it computes a multivariate Cox proportional hazard model using
#' all previously selected genes as predictor variables. Subsequently, it
#' calculates for each patient using this model their linear predictor,
#' which corresponds to the exponent of their proportionality constant.
#' @param genes_case_full_data Input matrix whose columns correspond to the patients
#' and rows to the genes, having selected only the columns belonging to disease
#' samples and having selected the rows corresponding to the genes selected
#' by the function \code{gene_selection_fun_filt}. The names of the rows must
#' be the names of the genes.
#' @param survival_time Numeric vector that includes time to the event
#' information after removal of NAs corresponding to healthy tissue samples
#' (case_tag == control_tag).
#' @param survival_event Numeric vector that indicates if relapse or death
#' have been produced (0 and 1s) after removal of NAs corresponding to healthy
#' tissue samples (case_tag == control_tag).
#' @return A numeric vector including the values produced by the function
#' for each sample in the dataset \code{genes_case_full_data}.
#' @export
#' @import survival
calculate_filter_function <- function(genes_case_full_data, survival_time,
                                      survival_event){
  obj_surv <-  survival::Surv(time = as.numeric(survival_time),
                              event = as.numeric(survival_event))
  obj_cox <- survival::coxph(obj_surv ~ .,
                             data = data.frame(t(genes_case_full_data)))

  filter_values <- obj_cox[["linear.predictors"]]
  names(filter_values) <- colnames(genes_case_full_data)

  return(filter_values)
}
