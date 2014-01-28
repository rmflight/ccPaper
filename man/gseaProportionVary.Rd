\name{gseaProportionVary}
\alias{gseaProportionVary}
\title{vary gsea proportion}
\usage{
  gseaProportionVary(inSize, sigProp, statistics)
}
\arguments{
  \item{inSize}{vector of sizes}

  \item{sigProp}{list of proportions, one for each sample}

  \item{statistics}{the statistics to use}
}
\value{
  list of lists for each \code{term}, with \emph{index} and
  \emph{stats}, \code{stats} contains a column for each
  sample, as well as a combined column.
}
\description{
  vary gsea proportion
}
\details{
  For each supplied inSize and sigProp, random indices of
  \code{inSize} are taken from the statistics. Based on
  \code{sigProp}, for each sample different numbers of
  indices take values from the upper 50 of ranked values
  and the lower 50 of ranked values. In this case each
  entry in \code{inSize} and denotes a different term, and
  they are treated independently. Each \emph{list} entry of
  \code{sigProp} is a sample, and each vector should be the
  same length as \code{inSize}. The \emph{list} combined
  method is based on averaging the statistics in this case.
}

