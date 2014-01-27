\name{goSample}
\alias{goSample}
\title{generate GO sample with genes}
\usage{
  goSample(go2gene, nGO, goLimits, nSample = 2,
    nGene = 1000, nNoise = 0, sharedNoiseFrac = 1)
}
\arguments{
  \item{go2gene}{\code{list} of go terms with associated
  genes}

  \item{nGO}{the number of go terms. see \emph{Details} for
  having different numbers of terms}

  \item{goLimits}{ranges of gene counts to restrict
  sampling. If multiple use a list}

  \item{nSample}{how many samples of genes do we want}

  \item{nGene}{the number of genes to sample that are
  annotated to the GO terms}

  \item{nNoise}{the number of noise genes}

  \item{sharedNoiseFrac}{how much of the noisy genes should
  be shared}
}
\value{
  list
}
\description{
  To test reproducibility of GO term sampling, we need to
  be able to sample GO terms and the genes annotated to
  them, and possibly noise genes as well.
}
\details{
  \code{nGO} and \code{goLimits} should be the same length
}

