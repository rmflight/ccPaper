\name{rankGenes}
\alias{rankGenes}
\title{rank genes using limma}
\usage{
  rankGenes(exprData, sampleStatus, doComps,
    dupStrategy = "minP", doAggregation = FALSE,
    aggregateIndex, adjust.method = "BH")
}
\arguments{
  \item{exprData}{log-expression values from an
  ExpressionSet}

  \item{pData}{}

  \item{sampleStatus}{character vector describing the
  samples}

  \item{doComps}{character vector of which comparisons to
  make}

  \item{dupStrategy}{how to resolve duplicates, default is
  to take the smallest p-value}

  \item{correction}{which multiple testing correction to
  apply to the results}
}
\value{
  list of data frames
}
\description{
  rank genes using limma
}

