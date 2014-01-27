\name{pvaluesMultiEnrich}
\alias{pvaluesMultiEnrich}
\title{get p-values from multiple enrichments}
\usage{
  pvaluesMultiEnrich(geneListNames, useTerms,
    hyperEnrichList, log = TRUE)
}
\arguments{
  \item{geneListNames}{the names to combine and take the
  minimum}

  \item{hyperEnrichList}{the set of hypergeometric
  enrichment results}

  \item{useTerms}{}

  \item{log}{return original or log-transformed values}
}
\value{
  \code{data.frame}, see details
}
\description{
  get p-values from multiple enrichments
}
\details{
  the \code{data.frame} returned has the p-values for each
  GO term from both the set based and list based results.
  The set-based value is the minimum value of all the set
  enrichments for that GO term, log transformed.
}

