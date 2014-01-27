\name{sampleTerms}
\alias{sampleTerms}
\title{samples genes from GO}
\usage{
  sampleTerms(terms, term2gene, nGenes = 2000, expRate = 4)
}
\arguments{
  \item{terms}{the set of GO terms to generate a sample
  for}

  \item{term2gene}{a list of the mapping of terms to genes}

  \item{nGenes}{the total number genes to sample from all
  genes annotated to the GO terms}
}
\value{
  list, with the full set of sampled genes in one entry,
  and a term2gene list with only those genes that are in
  the sample
}
\description{
  Generate random samples of genes from GO in such a way
  that GO terms with few annotations will have
  proportionately more genes than GO terms with lots of
  genes. This is accomplished by iterating over each GO
  term and adding sampled genes from the pool of available
  genes. The fraction of genes to add is based on a value
  from an exponential distribution.
}

