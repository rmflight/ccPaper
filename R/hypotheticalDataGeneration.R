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

#' calculate fraction of genes annotated to GO terms
#' 
#' For the supplied list of GO 2 gene annotations, and lists of genes, calculate the what fraction of a GO term is in the gene list.
#' 
#' @details returns a \code{data.frame} with the total count, and fraction for each of the gene lists supplied, suitable for plotting in \code{ggplot2}.
#' @param go2gene list of GO terms and genes they annotate
#' @param geneList list of character vectors of genes
#' @export
#' @return data.frame
calcFraction <- function(go2gene, geneList){
  nGO <- length(go2gene)
  
  goSizes <- sapply(go2gene, length)
  goNames <- names(go2gene)
  
  nList <- length(geneList)
  
  goFrac <- lapply(geneList, function(inGenes){
    sapply(go2gene, function(inGO){
      length(intersect(inGenes, inGO)) / length(inGO)
    })
  })
  goFrac <- stack(goFrac)
  names(goFrac) <- c("frac", "genelist")
  goFrac$goid <- rep(goNames, nList)
  goFrac$size <- rep(goSizes, nList)
  return(goFrac)
}

#' hyperGOMultiEnrichment
#' 
#' hypergeometric (Fisher's Exact Test) enrichment calculations for multiple lists of genes, with common background, as well as the intersection of the gene lists.
#' 
#' @param geneList \code{list} where each entry is a gene list
#' @param universe \code{character} vector of the background (assumed to be common across all)
#' @param ontology which GO ontology to use (default is "BP")
#' @param annotation what is the source of annotation (default is "org.Hs.eg.db")
#' @importClassesFrom categoryComparePaperRev GOHyperGParamsCC
#' @importFrom categoryComparePaperRev hyperGTestCC
#' @return \code{list}, see Details
#' @details the object returned will have an object for each sample list, as well as the one generated from the intersection
#' @export
hyperGOMultiEnrichment <- function(geneList, universe, ontology="BP", annotation="org.Hs.eg.db"){
  tmpIntersect <- geneList[[1]]
  
  nList <- length(geneList)
  for (iList in 1:nList){
    tmpIntersect <- intersect(tmpIntersect, geneList[[iList]])
  }
  
  geneList$intersect <- tmpIntersect
  
  outMulti <- lapply(geneList, function(inList){
    inHyper <- new("GOHyperGParamsCC", geneIds=inList, universeGeneIds=universe, ontology=ontology, fdr=0, annotation=annotation)
    outHyper <- hyperGTestCC(inHyper)
  })
}

#' get p-values from multiple enrichments
#' 
#' @param geneListNames the names to combine and take the minimum
#' @param hyperEnrichList the set of hypergeometric enrichment results
#' @param useTerms
#' @param log return original or log-transformed values
#' @export
#' @return \code{data.frame}, see details
#' @details the \code{data.frame} returned has the p-values for each GO term from both the set based and list based results. The set-based value is the minimum value of all the set enrichments for that GO term, log transformed.
pvaluesMultiEnrich <- function(geneListNames, useTerms, hyperEnrichList, log=TRUE){
  naVal <- 1
  if (log){
    naVal <- 0
  }
  
  # get the p-values
  allValues <- lapply(hyperEnrichList, function(inEnrich){
    tmpVal <- inEnrich@pvalues[useTerms]
    if (log){
      tmpVal <- -1 * log10(tmpVal)
    }
    names(tmpVal) <- useTerms
    tmpVal
  })
  
  setValues <- do.call(cbind, allValues[geneListNames])
  
  if (log){
    setValues <- apply(setValues, 1, min, na.rm=TRUE)
  } else {
    setValues <- apply(setValues, 1, max, na.rm=TRUE)
  }
  
  setValues[is.na(setValues)] <- naVal
  
  allValues$intersect[is.na(allValues$intersect)] <- naVal
  
  outValues <- data.frame(set = setValues, list = allValues$intersect)
  
  # get whether it was even measured
  presentValues <- lapply(hyperEnrichList, function(inEnrich){
    measuredGO <- names(inEnrich@pvalues)
    hasGO <- useTerms %in% measuredGO
    names(hasGO) <- useTerms
    hasGO
  })
  
  initPresent <- logical(length(useTerms))
  setPresent <- do.call(cbind, presentValues[geneListNames])
  setPresent <- apply(setPresent, 1, function(x){sum(x) > 0})
  
  outPresent <- data.frame(set = setPresent, list = presentValues$intersect)
  
  return(list(values = outValues, present = outPresent))
}

#' calculate difference and significance
#' 
#' @param pvalueData the \code{data.frame} of values
#' @param pCutoff the p-value cutoff to use for significance
#' @param log the values are logged, and p-value should be too
#' @return \code{data.frame}
#' @export
pvalueDiffSig <- function(pvalueData, pCutoff=0.05, log=TRUE){
  if (log){
    pCutoff <- -1 * log10(pCutoff)
  }
  
  pvalueData$diff <- pvalueData$set - pvalueData$list
  pvalueData$sigState <- "none"
  
  sigSets <- list(both = (pvalueData$set >= pCutoff) & (pvalueData$list >= pCutoff),
                  set = (pvalueData$set >= pCutoff) & (pvalueData$list < pCutoff),
                  list = (pvalueData$set < pCutoff) & (pvalueData$list >= pCutoff))
  
  invisible(lapply(names(sigSets), function(inSig){
    if (sum(sigSets[[inSig]]) > 0){
      pvalueData$sigState[sigSets[[inSig]]] <<- inSig
    }
  }))
  return(pvalueData)
}

#' trim go2gene and get fractions
#' 
#' for the cases of examining the "real" data, we need to take a GO 2 gene mapping, get the gene 2 GO mapping and find just those GO terms that actually have genes in the DE list, and finally calculate the fractions.
#' 
#' @param go2gene list of GO 2 gene mappings
#' @param deList the differentially expressed gene list as a \code{character} vector
#' @importFrom Biobase reverseSplit
#' @return \code{data.frame} of fractions and counts
#' @export
trimGOFrac <- function(go2gene, deList){
  stopifnot(is.list(go2gene))
  stopifnot(is.character(deList))
  
  go2gene <- lapply(go2gene, unique)
  
  gene2go <- reverseSplit(go2gene)
  gene2go <- gene2go[deList]
  keepGO <- unique(unlist(gene2go, use.names=F))
  
  go2gene <- go2gene[keepGO]
  
  goSize <- sapply(go2gene, length)
  goFrac <- sapply(go2gene, function(x){
    length(intersect(x, deList)) / length(x)
  })
  return(data.frame(size=goSize, frac=goFrac))
}


#' @name lung.RData
#' @title lung.RData
#' @docType data
#' @source ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18842/matrix/GSE18842_series_matrix.txt.gz
#' @details Downloaded, created an \code{ExpressionSet} object using \code{getGEO} on Dec 6, 2013.
NULL