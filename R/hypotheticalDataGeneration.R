#' returns random sample between limits
#' 
#' @param counts a named data.frame with the counts for each item (row)
#' @param limits the minimum and maximum limits to constrain the sample
#' @param nItem the number of items to return
#' @export
#' @return sampleNames character vector of row names sampled from within the limits provided
limitedRandomSample <- function(counts, limits, nItem=NULL){
  limitCounts <- counts[((counts[,1] >= limits[1]) & (counts[,1] < limits[2])), ,drop=F]
  if (!is.null(nItem)){
    sampleIndex <- sample(nrow(limitCounts), nItem)
  }
  else {
    sampleIndex <- seq(1, nrow(limitCounts))
  }
  return(rownames(limitCounts)[sampleIndex])
}

#' samples genes from GO
#' 
#' Generate random samples of genes from GO in such a way that GO terms with few annotations will have proportionately more genes than GO terms with lots of genes. This is accomplished by iterating over each GO term and adding sampled genes from the pool of available genes. The fraction of genes to add is based on a value from an exponential distribution.
#' 
#' @param terms the set of GO terms to generate a sample for
#' @param term2gene a list of the mapping of terms to genes
#' @param nGenes the total number genes to sample from all genes annotated to the GO terms
#' @export
#' @importFrom Biobase reverseSplit
#' @return list, with the full set of sampled genes in one entry, and a term2gene list with only those genes that are in the sample
sampleTerms <- function(terms, term2gene, nGenes=2000, expRate=4){
  term2gene <- term2gene[terms]
  gene2term <- reverseSplit(term2gene)
  
  allGenes <- names(gene2term)
  if (length(allGenes) < nGenes){
    nGenes <- length(allGenes)
  }
  sampleGenes <- character(nGenes)
  geneIndex <- 1
  haveSampled <- 0
  
  sampleTerms <- names(term2gene)
  sampleSizes <- sapply(term2gene, length)
  iTerm <- 1
  
  while((haveSampled <= nGenes) & (length(sampleTerms) > 0)){
    # do everything for one single term
    useTerm <- sampleTerms[iTerm]
    instanceSize <- ceiling(rexp(1, rate=expRate) * sampleSizes[useTerm])
    if(instanceSize > sampleSizes[useTerm]){
      instanceSize <- sampleSizes[useTerm]
    }
    instanceGenes <- sample(term2gene[[useTerm]], instanceSize)
    stopIndex <- geneIndex + instanceSize - 1
    sampleGenes[seq(geneIndex, stopIndex)] <- instanceGenes
    
    geneIndex <- stopIndex + 1 # set up for the next go around
    haveSampled <- haveSampled + instanceSize
    
    # now make any corrections needed to other terms and genes
    otherTerms <- unique(unlist(gene2term[instanceGenes]))
    otherSetDiff <- lapply(term2gene[otherTerms], function(x){
      setdiff(x, instanceGenes)
    })
    term2gene[otherTerms] <- otherSetDiff
    
    sampleSizes <- sapply(term2gene, length)
    not0 <- sampleSizes != 0
    
    sampleSizes <- sampleSizes[not0]
    term2gene <- term2gene[not0]
    sampleTerms <- names(term2gene)
    
    gene2term <- reverseSplit(term2gene)
    
    iTerm <- iTerm + 1
    if (iTerm >= length(term2gene)){
      iTerm <- 1
    } 
  }
  return(sampleGenes)
}