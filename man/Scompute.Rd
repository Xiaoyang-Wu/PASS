% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PASS_package.R
\name{Scompute}
\alias{Scompute}
\title{Computing the similarity matrix}
\usage{
Scompute(g, Z_test, d)
}
\arguments{
\item{g}{Character. The criterion function used for measuring similarity. Take "RBF" or "COS" and "RBF" as default.}

\item{Z_test}{Matrix. Normalized covariate matrix of the unlabeled data.}

\item{d}{Integral. Dimensionality of the covariates.}
}
\value{
The similarity matrix of the unlabeled data set.
}
\description{
Computing the similarity matrix
}
