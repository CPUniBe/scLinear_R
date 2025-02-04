% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sclinear.R
\name{predict}
\alias{predict}
\title{Predicts ADT levels based on the gene expression levels.}
\usage{
predict(predictor, gexp, layer = "counts", normalize_gex = TRUE)
}
\arguments{
\item{predictor}{predictor trained with 'fit_predictor' function}

\item{gexp}{Seurat Object containing gene expression data in the layer specified by argument layer}

\item{layer}{From which layer of the Seurat Object should the gene expression data be extracted}

\item{normalize_gex}{Should gene expression levels be normalized across cells (see function gexp_normalize for details)}
}
\value{
Predicted ADT levels. If ADT levels have been normalized during predictor training (default behaviour), the predicted ADT levels should also be considered 'normalized' in the same way.
}
\description{
Predicts ADT levels based on the gene expression levels.
}
