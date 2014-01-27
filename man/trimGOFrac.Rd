\name{trimGOFrac}
\alias{trimGOFrac}
\title{trim go2gene and get fractions}
\usage{
  trimGOFrac(go2gene, deList)
}
\arguments{
  \item{go2gene}{list of GO 2 gene mappings}

  \item{deList}{the differentially expressed gene list as a
  \code{character} vector}
}
\value{
  \code{data.frame} of fractions and counts
}
\description{
  for the cases of examining the "real" data, we need to
  take a GO 2 gene mapping, get the gene 2 GO mapping and
  find just those GO terms that actually have genes in the
  DE list, and finally calculate the fractions.
}

