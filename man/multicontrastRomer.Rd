\name{multicontrastRomer}
\alias{multicontrastRomer}
\title{multicontrast romer}
\usage{
  multicontrastRomer(index, y, design, contrastMatrix,
    nrot = 9999, do.parallel = TRUE)
}
\arguments{
  \item{index}{list of indices specifying how y maps to
  gene sets}

  \item{y}{numeric matrix of log-expression values}

  \item{design}{design matrix}

  \item{contrastMatrix}{the contrast.matrix, each column
  will be used}

  \item{nrot}{number of rotations}
}
\value{
  list of results for each contrast performed
}
\description{
  given a contrast.matrix with multiple columns, do
  \code{romer} for those
}
\seealso{
  romer
}

