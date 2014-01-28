\name{gseaProportion}
\alias{gseaProportion}
\title{generates proportion of genes in top 50}
\usage{
  gseaProportion(inSize, sigProp, statistics, nSample = 2)
}
\arguments{
  \item{inSize}{vector of sizes to sample}

  \item{sigProp}{vector of significant proportion, i.e.
  fraction above 50 rank}

  \item{statistics}{the set of statistics to use (will be
  ranked by function)}

  \item{nSample}{number of samples to use}
}
\value{
  list of lists for each \code{term}, with \emph{index} and
  \emph{stats}, \code{stats} contains a column for each
  sample, as well as a combined column.
}
\description{
  generates proportion of genes in top 50
}
\details{
  For each supplied inSize and sigProp, random indices of
  \code{inSize} are taken from the statistics. Based on
  \code{sigProp}, for each sample different numbers of
  indices take values from the upper 50 of ranked values
  and the lower 50 of ranked values. In this case each
  entry in \code{inSize} and \code{sigProp} denotes a
  different term, and they are treated independently. The
  \emph{list} combined method is based on averaging the
  statistics in this case.
}

