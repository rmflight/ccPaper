\name{calcFraction}
\alias{calcFraction}
\title{calculate fraction of genes annotated to GO terms}
\usage{
  calcFraction(go2gene, geneList)
}
\arguments{
  \item{go2gene}{list of GO terms and genes they annotate}

  \item{geneList}{list of character vectors of genes}
}
\value{
  data.frame
}
\description{
  For the supplied list of GO 2 gene annotations, and lists
  of genes, calculate the what fraction of a GO term is in
  the gene list.
}
\details{
  returns a \code{data.frame} with the total count, and
  fraction for each of the gene lists supplied, suitable
  for plotting in \code{ggplot2}.
}

