\name{multiSampleGeneSetTest}
\alias{multiSampleGeneSetTest}
\title{run geneSetTest for multiple samples}
\usage{
  multiSampleGeneSetTest(samplePValues, genesets,
    alternative = "mixed", ranks.only = FALSE,
    transform2Log = TRUE)
}
\arguments{
  \item{samplePValues}{list of p-values for each sample}

  \item{genesets}{list of gene sets to test}

  \item{alternative}{the alternative hypothesis to test}

  \item{ranks.only}{use the ranks only}
}
\value{
  list with the gene set p-values for each sample, as well
  as the combined
}
\description{
  For a list of samples and list of gene sets, do testing
  of the gene sets for each sample, as well as a
  meta-sample where there p-values for the samples are
  combined using Fisher's method
}
\seealso{
  geneSetTest
}

