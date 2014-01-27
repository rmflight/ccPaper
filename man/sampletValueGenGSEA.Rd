\name{sampletValueGenGSEA}
\alias{sampletValueGenGSEA}
\title{change set of t-values for GSEA}
\usage{
  sampletValueGenGSEA(sampleList, universeGenes, tValues)
}
\arguments{
  \item{sampleList}{the sample genes we will use}

  \item{universeGenes}{all the genes to use}

  \item{tValues}{a distribution of t-values}
}
\value{
  list of new t-values that should put genes from
  sampleList at the top
}
\description{
  Given a list of sample genes, the gene universe, and a
  set of t-values
}

