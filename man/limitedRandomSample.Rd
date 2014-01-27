\name{limitedRandomSample}
\alias{limitedRandomSample}
\title{returns random sample between limits}
\usage{
  limitedRandomSample(counts, limits, nItem = NULL)
}
\arguments{
  \item{counts}{a named data.frame with the counts for each
  item (row)}

  \item{limits}{the minimum and maximum limits to constrain
  the sample}

  \item{nItem}{the number of items to return}
}
\value{
  sampleNames character vector of row names sampled from
  within the limits provided
}
\description{
  returns random sample between limits
}

