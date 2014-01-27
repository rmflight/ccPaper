\name{gseaDiffs}
\alias{gseaDiffs}
\title{results and diffs following genesettest}
\usage{
  gseaDiffs(gseaRes, contrast = "down",
    orgLists = c("s1", "s2"))
}
\arguments{
  \item{gseaRes}{a list of results from doing geneSetTests}

  \item{contrast}{which alternative hypothesis to extract
  (for fake data only one should be valid)}

  \item{orgLists}{the names of the original samples to
  extract, will be compared with "comb"}
}
\value{
  data.frame of -1*log of p-values, and difference
}
\description{
  results and diffs following genesettest
}

