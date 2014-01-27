\name{hyperGOMultiEnrichment}
\alias{hyperGOMultiEnrichment}
\title{hyperGOMultiEnrichment}
\usage{
  hyperGOMultiEnrichment(geneList, universe,
    ontology = "BP", annotation = "org.Hs.eg.db")
}
\arguments{
  \item{geneList}{\code{list} where each entry is a gene
  list}

  \item{universe}{\code{character} vector of the background
  (assumed to be common across all)}

  \item{ontology}{which GO ontology to use (default is
  "BP")}

  \item{annotation}{what is the source of annotation
  (default is "org.Hs.eg.db")}
}
\value{
  \code{list}, see Details
}
\description{
  hypergeometric (Fisher's Exact Test) enrichment
  calculations for multiple lists of genes, with common
  background, as well as the intersection of the gene
  lists.
}
\details{
  the object returned will have an object for each sample
  list, as well as the one generated from the intersection
}

