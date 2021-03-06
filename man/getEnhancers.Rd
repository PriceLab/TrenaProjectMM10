% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TrenaProjectMM10.R
\name{getEnhancers,TrenaProjectMM10-method}
\alias{getEnhancers,TrenaProjectMM10-method}
\alias{getEnhancers}
\title{Get all the enhancer regions for the gene}
\usage{
\S4method{getEnhancers}{TrenaProjectMM10}(obj, targetGene = NA_character_, tissues = "all", maxSize = 10000)
}
\arguments{
\item{obj}{An object of class TrenaProjectMM10}

\item{targetGene}{default NA, in which case the current object's targetGene is used.}
}
\description{
Get all the enhancer regions for the gene
}
\seealso{
setTargetGene
}
