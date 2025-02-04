% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sclinear.R
\name{evaluate_predictor}
\alias{evaluate_predictor}
\title{Reports the mean RMSE of all ADT-specific regression models as well as the Pearson-/ & Spearman correlation coefficients
between predicted and measured ADT for each model separately
Unlike other predictor related functions, layers}
\usage{
evaluate_predictor(
  predictor,
  gexp_test,
  adt_test,
  gexp_layer = "counts",
  adt_layer = "counts",
  normalize_gex = TRUE,
  normalize_adt = TRUE,
  margin = 2
)
}
\arguments{
\item{predictor}{ADT-predictor (of form as produced by fit_predictor function) whose performance should be evaluated}

\item{gexp_test}{Seurat object containing gene expression data to use for ADT level prediction}

\item{adt_test}{Measured ADT levels to which the predicted ADT levels are to be compared}

\item{gexp_layer}{From which layer of the Seurat Object should the gene expression data be extracted}

\item{adt_layer}{From which layer of the Seurat Object should the gene expression data be extracted}

\item{normalize_gex}{Should gene expression levels be normalized across cells (see function gexp_normalize for details)}

\item{normalize_adt}{Should ADT levels be normalized (using the NormalizeData function from the Seurat package with method CLR) --> this needs
to be the same as the setting used during training of the prediction model.}

\item{margin}{Margin to apply CLR normalization over for ADT levels}
}
\value{
List of metrics for model performance (Mean RMSE of all models and ADT-specific correlation coefficients between predicted and measured ADT Values)
}
\description{
Reports the mean RMSE of all ADT-specific regression models as well as the Pearson-/ & Spearman correlation coefficients
between predicted and measured ADT for each model separately
Unlike other predictor related functions, layers
}
