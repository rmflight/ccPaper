
## ----nGO-----------------------------------------------------------------
nGO <- 20


## ----goAnnotations-------------------------------------------------------
library(org.Hs.eg.db)
library(GO.db)
library(ccPaperRev)
hsGO <- as.list(org.Hs.egGO2ALLEGS)

goOntology <- Ontology(names(hsGO))
goBP <- names(goOntology)[goOntology == "BP"]
hsGO <- hsGO[goBP]
hsGO <- lapply(hsGO, unique)
universeGenes <- unique(unlist(hsGO))


## ----goCountDist---------------------------------------------------------
library(ggplot2)
hsGO_count <- sapply(hsGO, length)
hsGO_count <- data.frame(count=hsGO_count)
ggplot(hsGO_count, aes(x=count)) + geom_bar(binwidth=10) + xlim(0, 2000) + ylim(0, 500)


## ----defineMin-----------------------------------------------------------
minGO <- 10


## ----defineGroupLimits---------------------------------------------------
grpLow <- c(10, 100)
grpMed <- c(250, 500)
grpHi <- c(500, 1500)


## ----sampleGroups--------------------------------------------------------
set.seed(271113) # for reproducibility
nGO <- 20

GO_low <- limitedRandomSample(hsGO_count, grpLow, nGO)
GO_med <- limitedRandomSample(hsGO_count, grpMed, nGO)
GO_hi <- limitedRandomSample(hsGO_count, grpHi, nGO)


## ----generateSamples-----------------------------------------------------
useGO <- c(GO_low, GO_med, GO_hi)
goList <- list(low=GO_low, med=GO_med, hi=GO_hi)
nGene <- 1000
sample1_org <- sampleTerms(useGO, hsGO, 1000, 4)[1:nGene]
sample2_org <- sampleTerms(useGO, hsGO, 1000, 4)[1:nGene]


## ----noiseGenes----------------------------------------------------------
nNoise <- 500
gene2HsGO <- reverseSplit(hsGO)
not_useGO <- sapply(gene2HsGO, function(x){
  sum(x %in% useGO) == 0
})
noiseGenes <- names(not_useGO)[not_useGO]
noiseGenes <- sample(noiseGenes, nNoise)

# check that we did this right, the fraction should not change after adding noise genes
useGO_frac <- sapply(hsGO[useGO], function(x){
  length(intersect(x, sample1_org)) / length(x)
})
sample1 <- c(sample1_org, noiseGenes)

useGO_fracNoise <- sapply(hsGO[useGO], function(x){
  length(intersect(x, sample1)) / length(x)
})
sample2 <- c(sample2_org, noiseGenes)
plot(useGO_frac, useGO_fracNoise)


## ----checkAssumptions----------------------------------------------------
# check assumptions
useGO_counts <- sapply(hsGO[useGO], length)
useGO_frac <- sapply(hsGO[useGO], function(x){
  length(intersect(x, sample1)) / length(x)
})
checkGO <- data.frame(counts=useGO_counts, frac=useGO_frac, sample="1")
useGO_frac2 <- sapply(hsGO[useGO], function(x){
  length(intersect(x, sample2)) / length(x)
})
checkGO <- rbind(checkGO, data.frame(counts=useGO_counts, frac=useGO_frac2, sample="2"))
ggplot(checkGO, aes(x=counts, y=frac, color=sample)) + geom_point()


## ----overlapDist---------------------------------------------------------
overlapDist <- sapply(useGO, function(x){
  baseAnn <- hsGO[[x]]
  t1 <- intersect(baseAnn, sample1)
  t2 <- intersect(baseAnn, sample2)
  c(length(baseAnn), length(t1), length(t2), length(intersect(t1,t2)))
})
overlapDist <- t(overlapDist)
overlapFrac <- c(overlapDist[,2] / overlapDist[,1], overlapDist[,3]/overlapDist[,1], overlapDist[,4]/overlapDist[,1])
overlapFrac <- data.frame(count=rep(overlapDist[,1], 3), 
                          frac=overlapFrac, 
                          id=c(rep("s1", 60), rep("s2", 60), rep("s1.s2", 60)))
overlapFrac
ggplot(overlapFrac, aes(x=count, y=frac, color=id)) + geom_point()


## ----ccInit--------------------------------------------------------------
s1GO <- new("GOHyperGParamsCC", geneIds=sample1, universeGeneIds=universeGenes, ontology="BP", fdr=0, annotation="org.Hs.eg.db")
s1_enrich <- hyperGTestCC(s1GO)

s2GO <- new("GOHyperGParamsCC", geneIds=sample2, universeGeneIds=universeGenes, ontology="BP", fdr=0, annotation="org.Hs.eg.db")
s2_enrich <- hyperGTestCC(s2GO)


## ----pvalueCut-----------------------------------------------------------
pval <- 0.05
minP <- -1 * log10(pval)


## ----compRes-------------------------------------------------------------
s1_goP <- -1 * log10((s1_enrich@pvalues)[useGO])
s2_goP <- -1 * log10((s2_enrich@pvalues)[useGO])

s_minP <- apply(cbind(s1_goP, s2_goP), 1, min)
sum(s_minP >= minP)


## ----intersectMethod-----------------------------------------------------
comGO <- new("GOHyperGParamsCC", geneIds=intersect(sample1, sample2), universeGeneIds=universeGenes, ontology="BP", fdr=0, annotation="org.Hs.eg.db")
com_enrich <- hyperGTestCC(comGO)
com_goP <- -1 * log10((com_enrich@pvalues)[useGO])


## ----numericComparison---------------------------------------------------
sum(s_minP >= minP)
sum(com_goP >= minP)

invisible(lapply(goList, function(x){
  print(sum(s_minP[x] >= minP))
  print(sum(com_goP[x] >= minP))
}))


## ----plotHistogram-------------------------------------------------------
s_count <- sapply(goList, function(x){
  sum(s_minP[x] >= minP)
})
com_count <- sapply(goList, function(x){
  sum(com_goP[x] >= minP)
})

countData <- data.frame(count=c(s_count, com_count), grp=rep(names(goList), 2), type=rep(c("s", "com"), each=3))
countData$grp <- factor(countData$grp, levels=c("low", "med", "hi"), ordered=TRUE)
ggplot(countData, aes(x=grp, y=count, fill=type)) + geom_bar(stat="identity", position="dodge") + ylim(0, 20)


## ----plotPDiffs----------------------------------------------------------
pDiff <- data.frame(pdiff = s_minP - com_goP,
                    annFrac = useGO_frac,
                    grp = rep(names(goList), each=20))
pDiff$grp <- factor(pDiff$grp, levels=c("low", "med", "hi"), ordered=T)

ggplot(pDiff, aes(x=annFrac, y=pdiff)) + geom_point()


## ----plotPDiffsGroup-----------------------------------------------------
ggplot(pDiff, aes(x=annFrac, y=pdiff)) + geom_point() + facet_grid(. ~ grp, scales="free_x")


## ----plotPDiffsMethod----------------------------------------------------
pDiff$sigGroup <- "none"
sigBoth <- names(s_minP)[(s_minP >= 1.3) & (com_goP >= 1.3)]
pDiff[sigBoth, "sigGroup"] <- "both"
sigS <- names(s_minP)[(s_minP >= 1.3) & (com_goP < 1.3)]
pDiff[sigS, "sigGroup"] <- "s"

ggplot(pDiff, aes(x=annFrac, y=pdiff, color=sigGroup)) + geom_point() + facet_grid(. ~ grp, scales="free_x")


