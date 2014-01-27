\name{sameGOStats}
\alias{sameGOStats}
\title{calculate statistics for all the same GO terms across multiple samples}
\usage{
  sameGOStats(inList)
}
\arguments{
  \item{inList}{the list of results to work with}
}
\value{
  data.frame
}
\description{
  For the case where we have calculated p-values and
  differences for random samples of genes, we need to
  calculate the mean and standard deviation for each GO
  term.
}
\details{
  expects each entry in \code{inList} to be a data.frame,
  and takes the column "diff" and does the statistics on
  it.
}

