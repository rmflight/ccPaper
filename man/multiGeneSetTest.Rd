\name{multiGeneSetTest}
\alias{multiGeneSetTest}
\title{run multiple geneSetTest's}
\usage{
  multiGeneSetTest(statistics, genesets,
    alternative = "down", ranks.only = FALSE)
}
\arguments{
  \item{statistics}{named vector of statistics}

  \item{genesets}{list of the gene sets to test}

  \item{alternative}{the alternative hypothesis to test}

  \item{ranks.only}{use the ranks of the entries only}
}
\value{
  set of p-values for the genesets
}
\description{
  run multiple geneSetTest's
}
\seealso{
  geneSetTest
}

