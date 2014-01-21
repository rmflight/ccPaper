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

#' calculate statistics for all the same GO terms across multiple samples
#' 
#' For the case where we have calculated p-values and differences for random samples of genes, we need to calculate the mean and standard deviation for each GO term.
#' 
#' @param inList the list of results to work with
#' @details expects each entry in \code{inList} to be a data.frame, and takes the column "diff" and does the statistics on it.
#' @return data.frame
#' @export
sameGOStats <- function(inList){
  getDiff <- function(inVar){
    return(inVar$diff)
  }
  
  allVals <- lapply(inList, getDiff)
  allVals <- do.call(cbind, allVals)
  meanVals <- apply(allVals, 1, mean, na.rm=T)
  sdVals <- apply(allVals, 1, sd, na.rm=T)
  minCI <- meanVals - 1.96 * (sdVals / sqrt(ncol(allVals)))
  maxCI <- meanVals + 1.96 * (sdVals / sqrt(ncol(allVals)))
  outVals <- data.frame(mean = meanVals, min = minCI, max = maxCI)
  return(outVals)
}


#' calculate difference between noise and clean differences for 100 gene samples
#' 
#' @param inList the list of results to work with
#' @return data.frame
#' @export
diffGOStats <- function(inList){
  getDiff <- function(inVar){
    diffDiff <- inVar$noise$diff - inVar$clean$diff
    return(diffDiff)
  }
  
  allVals <- lapply(inList, getDiff)
  allVals <- do.call(cbind, allVals)
  meanVals <- apply(allVals, 1, mean, na.rm=T)
  sdVals <- apply(allVals, 1, sd, na.rm=T)
  minCI <- meanVals - 1.96 * (sdVals / sqrt(ncol(allVals)))
  maxCI <- meanVals + 1.96 * (sdVals / sqrt(ncol(allVals)))
  outVals <- data.frame(mean = meanVals, min = minCI, max = maxCI)
  return(outVals)
}

#' generate GO sample with genes
#' 
#' To test reproducibility of GO term sampling, we need to be able to sample GO terms and the genes annotated to them, and possibly noise genes as well.
#' 
#' @param go2gene \code{list} of go terms with associated genes
#' @param nGO the number of go terms. see \emph{Details} for having different numbers of terms
#' @param goLimits ranges of gene counts to restrict sampling. If multiple use a list
#' @param nSample how many samples of genes do we want
#' @param nGene the number of genes to sample that are annotated to the GO terms
#' @param nNoise the number of noise genes
#' @param sharedNoiseFrac how much of the noisy genes should be shared
#' @export
#' @return list
#' @details \code{nGO} and \code{goLimits} should be the same length
goSample <- function(go2gene, nGO, goLimits, nSample=2, nGene=1000, nNoise=0, sharedNoiseFrac=1){
  stopifnot(names(nGO) == names(goLimits))
  goCount <- data.frame(count=sapply(go2gene, length))
  
  goSamples <- unlist(lapply(names(nGO), function(inName){
    limitedRandomSample(goCount, goLimits[[inName]], nGO[inName])
  }))
  
  goClass <- unlist(lapply(names(nGO), function(inName){rep(inName, nGO[inName])}))
  goSize <- goCount[goSamples, 'count']
  
  outSamples <- lapply(seq(1, nSample), function(iSample){
    sampleTerms(goSamples, go2gene, nGene)[1:nGene]
  })
  
  goSampleAnn <- sapply(goSamples, function(inGO){
    min(sapply(outSamples, function(inSample){
      length(intersect(go2gene[[inGO]], inSample))
    }))
  })
  
  goFrac <- goSampleAnn / goSize
  
  outGOData <- data.frame(id = as.character(goSamples), class = goClass, size = goSize, frac = goFrac)
  
  noiseSample <- vector("list", nSample)
  
  # and now the noise genes
  if (nNoise != 0){
    noiseGenes <- possibleNoise(go2gene, goSamples)
    
    nShare <- round(sharedNoiseFrac * nNoise)
    nInd <- nNoise - nShare
    
    shareNoise <- sample(noiseGenes, nShare)
    noiseGenes <- setdiff(noiseGenes, shareNoise)
    
    
    for (iSample in seq(1, nSample)){
      noiseSample[[iSample]] <- c(shareNoise, sample(noiseGenes, nInd))
      noiseGenes <- setdiff(noiseGenes, noiseSample[[iSample]])
    }
    
  }
  
  sampleNames <- paste("s", seq(1, nSample), sep="")
  names(outSamples) <- sampleNames
  names(noiseSample) <- sampleNames
  
  return(list(goData = outGOData, geneSample = outSamples, noiseSample = noiseSample))
}


#' possible noise genes to sample from
#' 
#' @param go2gene list of GO term to gene mapping
#' @param goSample character list of GO terms already sampled
#' @export
#' @return noise genes to sample from
possibleNoise <- function(go2gene, goSample){
  allGenes <- unique(unlist(go2gene))
  sampleGenes <- unique(unlist(go2gene[goSample]))
  posGenes <- allGenes[!(allGenes %in% sampleGenes)]
  return(posGenes)
}

#' Different values of noise
#' 
#' @param noiseGenes all possible noise genes we can sample from
#' @param nSamples how many samples do we need?
#' @param sizeNoise a numeric vector of how many noise genes we want to add (\code{integer})
#' @param fracShared a numeric vector of the fraction of shared noise genes to sweep over (\code{decimal})
#' @return \code{list} of \code{lists}
#' @details each entry in the list corresponds to a particular number of noise genes, and each list therein corresponds to the noise genes for each sample
#' @export
sweepNoiseSample <- function(noiseGenes, nSamples = 2, sizeNoise = seq(0, 1000, 10), fracShared = seq(0, 1, 0.01)){
  outSamples <- lapply(sizeNoise, function(inSize){
    sizeSamples <- lapply(fracShared, function(inFrac){
      sampleGenes <- vector('list', nSamples)
      
      tmpGenes <- noiseGenes # a gene variable we can modify as needed
      nShared <- round(inSize * inFrac)
      nIndependent <- inSize - nShared
      
      shareGenes <- sample(tmpGenes, nShared)
      
      tmpGenes <- tmpGenes[!(tmpGenes %in% shareGenes)]
      
      for (iSample in seq(1, nSamples)){
        #print(iSample)
        uniqGenes <- sample(tmpGenes, nIndependent)
        sampleGenes[[iSample]] <- c(shareGenes, uniqGenes)
        tmpGenes <- tmpGenes[!(tmpGenes %in% uniqGenes)]
      }
      return(sampleGenes)
    })
  })
  return(outSamples)
  
}

#' fishers method for combining p-values
#' 
#' @param x set of p-values to combine
#' @export
#' @return combined p-values
#' @details function from \url{http://mikelove.wordpress.com/2012/03/12/combining-p-values-fishers-method-sum-of-p-values-binomial/}
fishersMethod <- function(x){
  pchisq(-2 * sum(log(x)),df=2*length(x),lower=FALSE)
} 

#' change set of t-values for GSEA
#' 
#' Given a list of sample genes, the gene universe, and a set of t-values
#' 
#' @param sampleList the sample genes we will use
#' @param universeGenes all the genes to use
#' @param tValues a distribution of t-values
#' @return list of new t-values that should put genes from sampleList at the top
#' @export
sampletValueGenGSEA <- function(sampleList, universeGenes, tValues){
  tValues <- sort(tValues)
  nUniverse <- length(universeGenes)
  outT <- numeric(nUniverse)
  names(outT) <- universeGenes
  newTSample <- lapply(sampleList, function(inList){
    outTSample <- outT
    nGene <- length(inList)
    
    inList <- sample(inList, nGene)
    useT <- sample(tValues[1:nGene], nGene)
    
    outTSample[inList] <- useT
    
    otherT <- sample(tValues[seq(nGene+1, nUniverse)], nUniverse-nGene)
    otherGenes <- universeGenes[!(universeGenes %in% inList)]
    outTSample[otherGenes] <- otherT
    return(outTSample)
  })
  
  combSample <- do.call(cbind, newTSample)
  combSample <- rowMeans(combSample)
  newTSample$comb <- combSample
  return(newTSample)
}

#' run multiple geneSetTest's
#' 
#' @param statistics named vector of statistics
#' @param genesets list of the gene sets to test
#' @param alternative the alternative hypothesis to test
#' @param ranks.only use the ranks of the entries only 
#' @seealso geneSetTest
#' @importFrom limma geneSetTest
#' @export
#' @return set of p-values for the genesets
multiGeneSetTest <- function(statistics, genesets, alternative="down", ranks.only=FALSE){
  genesetStatistics <- sapply(genesets, function(inSet){
    index <- which(names(statistics) %in% inSet)
    geneSetTest(index, statistics, alternative=alternative, ranks.only=ranks.only)
  })
  return(genesetStatistics)
}

#' run geneSetTest for multiple samples
#' 
#' For a list of samples and list of gene sets, do testing of the gene sets for each sample, as well
#' as a meta-sample where there p-values for the samples are combined using Fisher's method
#' 
#' @param samplePValues list of p-values for each sample
#' @param genesets list of gene sets to test
#' @param alternative the alternative hypothesis to test
#' @param ranks.only use the ranks only
#' @seealso geneSetTest
#' @return list with the gene set p-values for each sample, as well as the combined
#' @export
multiSampleGeneSetTest <- function(samplePValues, genesets, alternative="mixed", ranks.only=FALSE, transform2Log=TRUE){
  pvalMatrix <- do.call(cbind, samplePValues)
  combPValues <- apply(pvalMatrix, 1, fishersMethod)
  samplePValues$comb <- combPValues
  if (transform2Log){
    samplePValues <- lapply(samplePValues, function(inVal){
      -1 * log(inVal)
    })
  }
  sampleSetValues <- lapply(samplePValues, multiGeneSetTest, genesets, alternative, ranks.only)
  return(sampleSetValues)
}


#' multiple-multiple geneSetTests
#' 
#' The default \code{geneSetTest} in \code{limma} will only do a single index and a single test. This function allows one to pass in a set of indices, and do all three tests (up, down, mixed).
#' 
#' @param index list of indices generated by \code{symbols2indices}
#' @param statistics the statistics to use
#' @param ranks.only whether to rank only, default is FALSE
#' @param ... other parameters for \code{geneSetTest}
#' 
#' @details Note that the function will try to do all three of up, down and mixed, but if there is no sign (i.e. F-statistic instead of t-statistic) the function will only do a test for "mixed"
#' @export
#' @return data.frame of p-values for each alternative test
#' @importFrom limma geneSetTest
mmGeneSetTest <- function(index, statistics, ranks.only=FALSE, ...){
  alternatives <- c("up", "down", "mixed")
  
  allsamesign <- all(statistics >= 0) || all(statistics <= 0)
  if (allsamesign){
    alternatives <- "mixed"
  }
  
  byIndex <- lapply(index, function(inIndices){
    byAlternative <- sapply(alternatives, function(inAlternative){
      geneSetTest(inIndices, statistics, inAlternative, ranks.only=ranks.only)
    })
  })
  
  pvalues <- do.call(rbind, byIndex)
  return(pvalues)
}

#' collapse probes to genes
#' 
#' take the median value of a set of probes to the associated genes in each sample. Make sure that \code{collapseBy} is in the same order as the \code{exprData} rows, and any rows you want removed should have \code{NA}.
#' 
#' @param exprData an expressionSet
#' @param collapseBy character vector of associated gene ids
#' @return exprData a matrix of median intensities
#' @export
collapseProbes <- function(exprData, collapseBy){
  naLoc <- is.na(collapseBy)
  exprData <- exprData[!naLoc,]
  collapseBy <- collapseBy[!naLoc]
  
  doMedian <- function(inValues){
    if (nrow(inValues) > 1){
      outValues <- apply(inValues, 2, median)
    } else {
      outValues <- inValues
    }
    return(outValues)
  }
  
  collapseData <- by(exprData, collapseBy, doMedian)
  collapseData <- do.call(rbind, collapseData)
  return(collapseData)
}

#' rank genes using limma
#'
#' @param exprData log-expression values from an ExpressionSet
#' @param pData 
#' @param sampleStatus character vector describing the samples
#' @param doComps character vector of which comparisons to make
#' @param dupStrategy how to resolve duplicates, default is to take the smallest p-value
#' @param correction which multiple testing correction to apply to the results
#' @return list of data frames
#' @export
#' @importFrom limma makeContrasts lmFit contrasts.fit eBayes topTable
rankGenes <- function(exprData, sampleStatus, doComps, dupStrategy="minP", doAggregation=FALSE, aggregateIndex, adjust.method="BH"){
  f <- factor(sampleStatus)
  design <- model.matrix(~0 + f)
  colnames(design) <- levels(f)
  
  contrast.matrix <- makeContrasts(contrasts=doComps, levels=design)
  fit <- lmFit(exprData, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  
  outData <- lapply(doComps, function(inComp){
    compData <- topTable(fit2, coef=inComp, number=Inf, adjust.method=adjust.method, sort.by="none")
    compData$id <- rownames(compData)
    if (doAggregation){
      naIndex <- is.na(aggregateIndex)
      compData$aggregateBy <- aggregateIndex
      compData$orgID <- rownames(compData)
      compData <- compData[!naIndex,]
      compData <- by(compData, compData$aggregateBy, function(inData){
        if (dupStrategy == "minP"){
          inData[which.min(inData$P.Value),]
        }
        
      })
      compData <- do.call(rbind, compData)
    }
    return(compData)
  })
  return(outData)
}

#' get differentially expressed genes in a data.frame
#' 
#' @param rankedData data.frame of p-values and fold changes
#' @param id the column to use as a returned id
#' @param useP the column with the p-values
#' @param pcutoff maximum p-value to use
#' @param lfc column identifier for log-fold-change
#' @param lfcCutoff should a minumum fold-change be required
#' @param splitDir should significant results be split by up and down?
#' @export
#' @return list
#' @details if \code{splitDir=TRUE} (default), then the returned list will have \code{up}, \code{dn}, and \code{universe}. If \code{splitDir=FALSE}, then it will be \code{sig} and \code{universe}.
getDiffGenes <- function(rankedData, id="aggregateBy", useP="adj.P.Val", pcutoff=0.05, lfc="logFC", lfcCutoff=NA, splitDir=TRUE){
  sigEntry <- rankedData[, useP] <= pcutoff
  
  if (!(is.na(lfcCutoff))){
    sigEntry <- sigEntry & abs(rankedData[, lfc] >= lfcCutoff)
  }
  
  if (splitDir){
    upEntry <- rankedData[, lfc] > 0
    upID <- rankedData[(upEntry & sigEntry), id]
    dnEntry <- rankedData[, lfc] < 0
    dnID <- rankedData[(dnEntry & sigEntry), id]
    
    outData <- list(up=upID, dn=dnID, universe=rankedData[, id])
  } else {
    sigID <- rankedData[sigEntry, id]
    outData <- list(sig=sigID, universe=rankedData[, id])
  }
  
  return(outData)
  
}


#' results and diffs following genesettest
#' 
#' @param gseaRes a list of results from doing geneSetTests
#' @param contrast which alternative hypothesis to extract (for fake data only one should be valid)
#' @param orgLists the names of the original samples to extract, will be compared with "comb"
#' @export
#' @return data.frame of -1*log of p-values, and difference
gseaDiffs <- function(gseaRes, contrast="down", orgLists=c("s1", "s2")){
  orgRes <- do.call(cbind, lapply(gseaRes[orgLists], function(x){-1*log(x[,contrast])}))
  
  useNames <- rownames(orgRes)
  
  orgRes <- rowMin(orgRes)
  
  outRes <- data.frame(org=orgRes, comb=-1*log(gseaRes$comb[,contrast]))
  outRes$diff <- outRes$org - outRes$comb
  return(outRes)
}


#' multicontrast romer
#' 
#' given a contrast.matrix with multiple columns, do \code{romer} for those
#' 
#' @param index list of indices specifying how y maps to gene sets
#' @param y numeric matrix of log-expression values
#' @param design design matrix
#' @param contrastMatrix the contrast.matrix, each column will be used
#' @param nrot number of rotations
#' @seealso romer
#' @export
#' @return list of results for each contrast performed
#' @importFrom limma romer
multicontrastRomer <- function(index, y, design, contrastMatrix, nrot=9999){
  nContrast <- ncol(contrastMatrix)
  
  outData <- lapply(seq(1, nContrast), function(inCol){
    useContrast <- contrastMatrix[, inCol]
    romer(index, y, design, contrast=contrastMatrix[, inCol], nrot=nrot)
  })
}

#' generate entries within a distribution
#' 
#' @param inLimits the range to use
#' @param nSample how many entries within the range
#' 
#' @return set of sizes
#' @export
gseaSizeRange <- function(inLimits, nSample){
  round(runif(nSample, inLimits[1], inLimits[2]))
}

#' generates proportion of genes in top 50
#' 
#' @param inSize vector of sizes to sample
#' @param sigProp vector of significant proportion, i.e. fraction above 50 rank
#' @param statistics the set of statistics to use (will be ranked by function)
#' @param nSample number of samples to use
#' 
#' @details For each supplied inSize and sigProp, random indices of \code{inSize} are taken from the statistics. Based on \code{sigProp}, for each sample different numbers of indices take values from the upper 50 of ranked values and the lower 50 of ranked values. In this case each entry in \code{inSize} and \code{sigProp} denotes a different term, and they are treated independently. The \emph{list} combined method is based on averaging the statistics in this case.
#' 
#' @export
#' @return list of lists for each \code{term}, with \emph{index} and \emph{stats}, \code{stats} contains a column for each sample, as well as a combined column. 
gseaProportion <- function(inSize, sigProp, statistics, nSample=2){
  statistics <- sort(statistics)
  nPosGene <- length(statistics)
  grpLoc <- seq(1, nPosGene, round(nPosGene / 4))
  topStats <- statistics[seq(1, grpLoc[3])]
  botStats <- statistics[seq(grpLoc[3], nPosGene)]
  gseaStatistics <- lapply(seq(1, length(inSize)), function(inIndex){
    nGene <- inSize[inIndex]
    outIndex <- sample(nPosGene, nGene)
    useProp <- sigProp[inIndex]
    tmpStat <- statistics
    
    nTop <- round(nGene*useProp)
    nBot <- nGene - nTop
    
    outStatistics <- lapply(seq(1, nSample), function(inSample){
      topIndex <- sample(outIndex, nTop)
      botIndex <- outIndex[!(outIndex %in% topIndex)]
      tmpT <- tmpStat
      tmpT[topIndex] <- sample(topStats, nTop)
      tmpT[botIndex] <- sample(botStats, nBot)
      return(tmpT)
    })
    
    outStatistics <- do.call(cbind,outStatistics)
    outStatistics <- cbind(outStatistics, rowMeans(outStatistics))
    colnames(outStatistics) <- c(paste("s", seq(1, nSample), sep=""), "comb")
    
    return(list(index=outIndex, stats=outStatistics))
  })
  return(gseaStatistics)
}



#' vary gsea proportion
#' 
#' @param inSize vector of sizes
#' @param sigProp list of proportions, one for each sample
#' @param statistics the statistics to use
#' 
#' @details For each supplied inSize and sigProp, random indices of \code{inSize} are taken from the statistics. Based on \code{sigProp}, for each sample different numbers of indices take values from the upper 50 of ranked values and the lower 50 of ranked values. In this case each entry in \code{inSize} and denotes a different term, and they are treated independently. Each \emph{list} entry of \code{sigProp} is a sample, and each vector should be the same length as \code{inSize}. The \emph{list} combined method is based on averaging the statistics in this case.
#' 
#' @export
#' @return list of lists for each \code{term}, with \emph{index} and \emph{stats}, \code{stats} contains a column for each sample, as well as a combined column. 
gseaProportionVary <- function(inSize, sigProp, statistics){
  statistics <- sort(statistics)
  nPosGene <- length(statistics)
  nSample <- length(sigProp)
  grpLoc <- seq(1, nPosGene, round(nPosGene / 4))
  topStats <- statistics[seq(1, grpLoc[4])]
  botStats <- statistics[seq(grpLoc[4], nPosGene)]
  gseaStatistics <- lapply(seq(1, length(inSize)), function(inIndex){
    
    nGene <- inSize[inIndex]
    outIndex <- sample(nPosGene, nGene)
    tmpStat <- statistics
    
    outStatistics <- lapply(seq(1, nSample), function(inSample){
      # get significance proportion for each sample individually
      nTop <- round(nGene * sigProp[[inSample]][inIndex])
      nBot <- nGene - nTop
      topIndex <- sample(outIndex, nTop)
      botIndex <- outIndex[!(outIndex %in% topIndex)]
      tmpT <- tmpStat
      tmpT[topIndex] <- sample(topStats, nTop)
      tmpT[botIndex] <- sample(botStats, nBot)
      return(tmpT)
    })
    
    outStatistics <- do.call(cbind,outStatistics)
    outStatistics <- cbind(outStatistics, rowMeans(outStatistics))
    colnames(outStatistics) <- c(paste("s", seq(1, nSample), sep=""), "comb")
    
    return(list(index=outIndex, stats=outStatistics))
  })
  return(gseaStatistics)
}

#' apply regression model to pvalues
#' 
#' To transform p-values to a signed range, we will apply the results of a regression model to the p-values
#' 
#' @param inModel result from \code{lm}
#' @param pvalues set of pvalues to transform
#' @export 
#' @return transformed values
p2signed <- function(inModel, pvalues){
  outSigned <- inModel$coefficients[1] + inModel$coefficients[2]*pvalues
}

#' gseaProportion, combined using Fisher's method
#' 
#' @param inSize vector of sizes to sample
#' @param sigProp vector of significant proportion, i.e. fraction above 50 rank
#' @param statistics the set of statistics to use (will be ranked by function)
#' @param nSample number of samples to use
#' 
#' @details For each supplied inSize and sigProp, random indices of \code{inSize} are taken from the statistics. Based on \code{sigProp}, for each sample different numbers of indices take values from the upper 50 of ranked values and the lower 50 of ranked values. In this case each entry in \code{inSize} and \code{sigProp} denotes a different term, and they are treated independently. The \emph{list} combined method is based on averaging the statistics in this case.
#' 
#' @export
#' @return list of lists for each \code{term}, with \emph{index} and \emph{stats}, \code{stats} contains a column for each sample, as well as a combined column. 
gseaProportionFisher <- function(inSize, sigProp, statistics, nSample=2){
  statistics <- sort(statistics)
  nPosGene <- length(statistics)
  grpLoc <- seq(1, nPosGene, round(nPosGene / 4))
  topStats <- statistics[seq(1, grpLoc[3])]
  botStats <- statistics[seq(grpLoc[3]+1, nPosGene)]
  gseaStatistics <- lapply(seq(1, length(inSize)), function(inIndex){
    
    nGene <- inSize[inIndex]
    outIndex <- sample(nPosGene, nGene)
    tmpStat <- statistics
    
    outStatistics <- lapply(seq(1, nSample), function(inSample){
      # get significance proportion for each sample individually
      nTop <- round(nGene * sigProp[inIndex])
      nBot <- nGene - nTop
      topIndex <- sample(outIndex, nTop)
      botIndex <- outIndex[!(outIndex %in% topIndex)]
      tmpT <- tmpStat
      tmpT[topIndex] <- sample(topStats, nTop)
      tmpT[botIndex] <- sample(botStats, nBot)
      return(tmpT)
    })
    
    outStatistics <- do.call(cbind,outStatistics)
    combStat <- apply(outStatistics, 1, fishersMethod)
    outStatistics <- cbind(outStatistics, combStat)
    colnames(outStatistics) <- c(paste("s", seq(1, nSample), sep=""), "comb")
    
    return(list(index=outIndex, stats=outStatistics))
  })
  return(gseaStatistics)
}

#' @name lung.RData
#' @title lung.RData
#' @docType data
#' @source ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE18nnn/GSE18842/matrix/GSE18842_series_matrix.txt.gz
#' @details Downloaded, created an \code{ExpressionSet} object using \code{getGEO} on Dec 6, 2013.
NULL