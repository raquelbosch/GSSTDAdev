% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSSTDA_objects.R
\name{gene_selection_}
\alias{gene_selection_}
\title{Private gene_selection_}
\usage{
gene_selection_(
  full_data,
  survival_time,
  survival_event,
  control_tag_cases,
  gene_select_surv_type,
  num_gen_select_for_fun_filt,
  gene_select_mapper_metric,
  num_gen_select_for_mapper,
  matrix_disease_component = NULL
)
}
\arguments{
\item{full_data}{Input matrix whose columns correspond to the patients and
rows to the genes.}

\item{survival_time}{Numerical vector of the same length as the number of
columns of \code{full_data}. In addition, the patients must be in the same
order as in \code{full_data}. For the patients whose sample is pathological
should be indicated the time between the disease diagnosis and event
(death, relapse or other). If the event has not occurred, it should be
indicated the time until the end of follow-up. Patients whose sample is
from healthy tissue must have an NA value}

\item{survival_event}{Numerical vector of the same length as the number of
columns of \code{full_data}. Patients must be in the same order as in
\code{full_data}. For the the patients with pathological sample should
be indicated whether the event has occurred (1) or not (0). Only these
values are valid and healthy patients must have an NA value.}

\item{control_tag_cases}{Numeric vector with the indices of the columns of
\code{full_data} and/or \code{matrix_disease_component}
corresponding to the healthy sample patients.}

\item{gene_select_surv_type}{Option. Options on how to select the genes to be
used in the calculation of the values of the filter function (and in the gene
selection for Mapper if the option "sd_surv" has been chosen). Select the
"Abs" option, which means that the genes with the highest absolute value are
chosen, or the "Top_Bot" option, which means that half of the selected genes
are those with the highest value (positive value, i.e. worst survival
prognosis) and the other half are those with the lowest value (negative
value, i.e. best prognosis). "Top_Bot" default option.}

\item{num_gen_select_for_fun_filt}{Number of genes to be selected to be used
for the calculation of the values of the filter function.}

\item{gene_select_mapper_metric}{Gene selection criteria for Mapper. Choose as
selection criteria between:
\itemize{
\item "mad": those genes in the disease component matrix with the highest mean
absolute deviation will be selected
\item "sd": those genes with the highest standard deviation will be selected
\item “iqr": those genes with the highest interquartile range will be selected
\item "mean_sd": to choose genes with high mean and standard deviation simultaneously
\item "sd_surv": these option selects genes for mapper based on the product of standard
deviation of the genes in the disease component matrix plus one times
the Z score obtained by fitting a cox proportional hazard model to the level
of each gene.
}}

\item{num_gen_select_for_mapper}{Number of genes to be selected to be used
in mapper.}

\item{matrix_disease_component}{Optional, only necessary in case of gene
selection after DSGA has been performed. Matrix of the disease components
(the transformed \code{full_data} matrix from which the normal component has
been removed) from the \code{dsga_function}.}
}
\value{
A \code{gene_selection_object}. It contains:
\itemize{
\item the matrix with which the gene selection has been performed without NAN's
values (\code{data}). It is the \code{matrix_disease_component} in case it has been
performed from a \code{dsga_object} or \code{full_data} in the opposite case.
\item the \code{cox_all_matrix} (a matrix with the results of the application of
proportional hazard models: with the regression coefficients, the odds ratios,
the standard errors of each coefficient, the Z values (coef/se_coef) and
the p-values for each Z value)
\item a vector with the name of the selected genes for mapper
(\code{genes_selected_mapper}) and for the calculation of the values of the
filter function (\code{genes_selected_fun_filt}).
\item the matrix of disease components with only the rows of the selected genes
(\code{genes_disease_component})
\item and the vector of the values of the filter function of pathological samples
\code{filter_values}.
}
}
\description{
Private function to gene selection
}
\examples{
\donttest{
gene_select_surv_type <- "Top_Bot"
percent_gen_select <- 10
control_tag_cases <- which(case_tag == "NT")
gene_selection_obj <- gene_selection_(full_data, survival_time, survival_event,
control_tag_cases, gene_select_surv_type ="top_bot", num_gen_select_for_fun_filt = 200,
gene_select_mapper_metric="mad", num_gen_select_for_mapper = 1000)}
}
