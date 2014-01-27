\name{samplePValueGenGSEA}
\alias{samplePValueGenGSEA}
\title{change set of p-values for GSEA}
\usage{
  samplePValueGenGSEA(sampleList, initPValues, distPValues)
}
\arguments{
  \item{sampleList}{the sample genes we will use}

  \item{initPValues}{initial p-values for all the genes
  (named vector)}

  \item{distPValues}{a distribution of p-values to take
  samples from for the sample genes}
}
\value{
  list where for each sample the appropriate entries have
  modified p-values
}
\description{
  Given a list of sample genes, the full list of named
  p-values, and a mean and sd
}

