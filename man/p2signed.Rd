\name{p2signed}
\alias{p2signed}
\title{apply regression model to pvalues}
\usage{
  p2signed(inModel, pvalues)
}
\arguments{
  \item{inModel}{result from \code{lm}}

  \item{pvalues}{set of pvalues to transform}
}
\value{
  transformed values
}
\description{
  To transform p-values to a signed range, we will apply
  the results of a regression model to the p-values
}

