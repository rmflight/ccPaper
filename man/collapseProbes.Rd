\name{collapseProbes}
\alias{collapseProbes}
\title{collapse probes to genes}
\usage{
  collapseProbes(exprData, collapseBy)
}
\arguments{
  \item{exprData}{an expressionSet}

  \item{collapseBy}{character vector of associated gene
  ids}
}
\value{
  exprData a matrix of median intensities
}
\description{
  take the median value of a set of probes to the
  associated genes in each sample. Make sure that
  \code{collapseBy} is in the same order as the
  \code{exprData} rows, and any rows you want removed
  should have \code{NA}.
}

