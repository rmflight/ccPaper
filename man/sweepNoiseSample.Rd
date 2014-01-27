\name{sweepNoiseSample}
\alias{sweepNoiseSample}
\title{Different values of noise}
\usage{
  sweepNoiseSample(noiseGenes, nSamples = 2,
    sizeNoise = seq(0, 1000, 10),
    fracShared = seq(0, 1, 0.01))
}
\arguments{
  \item{noiseGenes}{all possible noise genes we can sample
  from}

  \item{nSamples}{how many samples do we need?}

  \item{sizeNoise}{a numeric vector of how many noise genes
  we want to add (\code{integer})}

  \item{fracShared}{a numeric vector of the fraction of
  shared noise genes to sweep over (\code{decimal})}
}
\value{
  \code{list} of \code{lists}
}
\description{
  Different values of noise
}
\details{
  each entry in the list corresponds to a particular number
  of noise genes, and each list therein corresponds to the
  noise genes for each sample
}

