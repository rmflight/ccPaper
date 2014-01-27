\name{getDiffGenes}
\alias{getDiffGenes}
\title{get differentially expressed genes in a data.frame}
\usage{
  getDiffGenes(rankedData, id = "aggregateBy",
    useP = "adj.P.Val", pcutoff = 0.05, lfc = "logFC",
    lfcCutoff = NA, splitDir = TRUE)
}
\arguments{
  \item{rankedData}{data.frame of p-values and fold
  changes}

  \item{id}{the column to use as a returned id}

  \item{useP}{the column with the p-values}

  \item{pcutoff}{maximum p-value to use}

  \item{lfc}{column identifier for log-fold-change}

  \item{lfcCutoff}{should a minumum fold-change be
  required}

  \item{splitDir}{should significant results be split by up
  and down?}
}
\value{
  list
}
\description{
  get differentially expressed genes in a data.frame
}
\details{
  if \code{splitDir=TRUE} (default), then the returned list
  will have \code{up}, \code{dn}, and \code{universe}. If
  \code{splitDir=FALSE}, then it will be \code{sig} and
  \code{universe}.
}

