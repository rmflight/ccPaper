
## ----customCSS, include=FALSE--------------------------------------------
cssFile <- system.file("extdata", "style.css", package="ccPaper")
options(markdown.HTML.stylesheet = cssFile)


## ----setupPlot-----------------------------------------------------------
library(ggplot2)
useTheme <- theme_bw() + theme(legend.key=element_blank(), legend.title = element_text(size=25), axis.title = element_text(size=25), axis.text = element_text(size=15), legend.text = element_text(size=20), strip.text.y = element_text(size=20), strip.text.x = element_text(size=20))


savePlot <- function(plotObj, plotFile, figPath="/mlab/data/rmflight/Documents/projects/work/ccPaper/savedPlots"){
  if (file.exists(figPath)){
    plotObj <- plotObj + ggtitle(NULL)
    png(filename=file.path(figPath, plotFile), width=854, height=628)
    print(plotObj)
    dev.off()
  }
}


## ----goAnnotations, message=FALSE----------------------------------------
library(org.Hs.eg.db)
library(GO.db)
library(ccPaper)
hsGO <- as.list(org.Hs.egGO2ALLEGS)

goOntology <- Ontology(names(hsGO))
goBP <- names(goOntology)[goOntology == "BP"]
hsGO <- hsGO[goBP]
hsGO <- lapply(hsGO, unique)
universeGenes <- unique(unlist(hsGO))


## ----goCountDist---------------------------------------------------------
hsGO_count <- sapply(hsGO, length)
hsGO_count <- data.frame(size=hsGO_count)
p <- ggplot(hsGO_count, aes(x=size)) + geom_bar(binwidth=10) + xlim(0, 2000) + ylim(0, 500) + useTheme + ggtitle("Number of genes annotated to GO terms") + xlab('Num. of Genes') + ylab('Count')
p

savePlot(p, "go2gene_distribution.png")


## ----defineGroupLimits---------------------------------------------------
grpLow <- c(10, 100)
grpMed <- c(250, 500)
grpHi <- c(500, 1500)


## ----sampleGroups--------------------------------------------------------
set.seed(271113) # for reproducibility
nGO <- c(50, 30, 20)
names(nGO) <- c("low", "med", "hi")

GO_low <- limitedRandomSample(hsGO_count, grpLow, nGO["low"])
GO_med <- limitedRandomSample(hsGO_count, grpMed, nGO["med"])
GO_hi <- limitedRandomSample(hsGO_count, grpHi, nGO["hi"])

useGO <- c(GO_low, GO_med, GO_hi)
useGOSizeClass <- c(rep("low", nGO["low"]), rep("med", nGO["med"]), rep("hi", nGO["hi"]))
names(useGOSizeClass) <- useGO


## ----deALL---------------------------------------------------------------
library(ALL)
library(limma)
data(ALL)

grpALL <- factor(strtrim(pData(ALL)$BT, 1))
designALL <- model.matrix(~0+grpALL)
colnames(designALL) <- c("B", "T")
fitALL <- lmFit(ALL, designALL)
contMatrixALL <- makeContrasts(BvT=B-T, levels=designALL)
fitALL2 <- contrasts.fit(fitALL, contMatrixALL)
fitALL2 <- eBayes(fitALL2)
deALL <- rownames(topTable(fitALL2, adjust="BH", p.value=0.05, number=Inf))


## ----goDistALL-----------------------------------------------------------
library(hgu95av2.db)
ALLgo <- as.list(hgu95av2GO2ALLPROBES)
ALL_sizeFrac <- trimGOFrac(ALLgo, deALL)
ggplot(ALL_sizeFrac, aes(x=size, y=frac)) + geom_point() + xlim(0, 1500) + useTheme + ggtitle("DE ALL fraction as a function of size") + xlab("Number of Genes Annotated to Term") + ylab("Fraction of Term Genes in DE Sample")


## ----deLung--------------------------------------------------------------
data(lung)
sampleStat <- strsplit2(pData(lung)$title, " ")

notNP <- grep("^NP", sampleStat[,2], invert=T)
lung <- lung[, notNP]
sampleStat <- sampleStat[notNP,]

sampleID <- factor(sampleStat[,2])
cancerID <- factor(sampleStat[,1])

lungDesign <- model.matrix(~sampleID+cancerID)
fitLung <- lmFit(lung, lungDesign)
fitLung <- eBayes(fitLung)
deLung <- rownames(topTable(fitLung, coef="cancerIDTumor", number=2000, adjust.method="BY", p.value=0.00001))


## ----goDistLung----------------------------------------------------------
library(hgu133plus2.db)
lung2go <- as.list(hgu133plus2GO2ALLPROBES)
lung_sizeFrac <- trimGOFrac(lung2go, deLung)
lung_sizeFrac$genelist <- "lung"
ggplot(lung_sizeFrac, aes(x=size, y=frac)) + geom_point() + xlim(0, 1500) + useTheme + ggtitle("DE Lung fraction as a function of size") + xlab("Num. of Genes") + ylab("Fraction")


## ----generateSamples-----------------------------------------------------
goList <- list(low=GO_low, med=GO_med, hi=GO_hi)
nGene <- 1000
sample1_org <- sampleTerms(useGO, hsGO, 1000, 4)[1:nGene]
sample2_org <- sampleTerms(useGO, hsGO, 1000, 4)[1:nGene]

#check our expectations
goFractions <- calcFraction(hsGO[useGO], list(sample1=sample1_org, sample2=sample2_org))
ggplot(goFractions, aes(x=size, y=frac, color=genelist)) + geom_point() + useTheme + ggtitle("Sample GO fractions") + xlab("Num. of Genes") + ylab("Fraction")


## ----compareFracs--------------------------------------------------------
sampleGOFracs <- calcFraction(hsGO, list(sample1=sample1_org, sample2=sample2_org))
ALL_sizeFrac$genelist <- "all"
ALL_tmp <- rbind(ALL_sizeFrac, lung_sizeFrac, sampleGOFracs[,c("frac", "genelist", "size")])
ggplot(ALL_tmp, aes(x=size, y=frac)) + geom_point() + xlim(0, 1500) + facet_grid(genelist ~ .) + useTheme + ggtitle("Sample GO fractions for all") + xlab("Num. of Genes") + ylab("Fraction")


## ----cleanupData---------------------------------------------------------
keepItems <- c("goFractions", "sample1_org", "sample2_org", "hsGO", "useGO", "useGOSizeClass", "GO_hi", "GO_low", "GO_med", "grpHi", "grpLow", "grpMed", "minGO", "nGO", "nGene", "universeGenes", "savePlot", "useTheme", "fitALL2")
allItems <- ls()
delItems <- allItems[!(allItems %in% keepItems)]
rm(list=delItems)


## ----runCalcs------------------------------------------------------------
samples_noNoise <- list(sample1 = sample1_org, sample2 = sample2_org)
go_noNoise <- hyperGOMultiEnrichment(samples_noNoise, universe=universeGenes)


## ----noNoisePvalues------------------------------------------------------
noNoise_results <- pvaluesMultiEnrich(c("sample1", "sample2"), useGO, go_noNoise)

noNoise_pvalues <- pvalueDiffSig(noNoise_results$values, pCutoff=0.05, log=TRUE)


## ----noNoiseAddInfo------------------------------------------------------
noNoise_pvalues$sizeClass <- useGOSizeClass
noNoise_pvalues$size <- goFractions[1:length(useGOSizeClass),4]
noNoise_pvalues$frac <- goFractions$frac[1:length(useGOSizeClass)]


## ----noNoisePlotStuff----------------------------------------------------
noNoise_pvalues$sizeClass <- factor(noNoise_pvalues$sizeClass, levels=c("low", "med", "hi"), ordered=TRUE)
p <- ggplot(noNoise_pvalues, aes(x=frac, y=diff)) + geom_point() + facet_grid(. ~ sizeClass, scales="free_x") + xlab("fraction") + useTheme + ggtitle("Sample - Combined Differences") + xlab("Fraction") + ylab("Diff. log p-values")
p
savePlot(p, "nonoise_sample_combineddiff.png")
ggplot(noNoise_pvalues, aes(x=sizeClass, y=diff)) + geom_boxplot() + geom_point() + ggtitle("Sample - Combined Differences") + xlab("Size Class") + ylab("Diff. log p-values")


## ----noNoiseBetterPVals--------------------------------------------------
tapply(noNoise_pvalues$diff, noNoise_pvalues$sizeClass, function(x){sum(x > 0)})


## ----noNoiseBetterMean---------------------------------------------------
tapply(noNoise_pvalues$diff, noNoise_pvalues$sizeClass, mean)


## ----noiseGenes----------------------------------------------------------
nNoise <- 500
noiseGenes <- possibleNoise(hsGO, useGO)
noiseGenes <- sample(noiseGenes, nNoise)

# check that we did this right, the fraction should not change after adding noise genes
sample1 <- c(sample1_org, noiseGenes)
sample2 <- c(sample2_org, noiseGenes)

samplesNoise <- list(sample1=sample1, sample2=sample2)

goFracNoise <- calcFraction(hsGO[useGO], samplesNoise)
plot(goFracNoise$frac, goFractions$frac)


## ----noisyEnrichment-----------------------------------------------------
go_noise <- hyperGOMultiEnrichment(samplesNoise, universeGenes)


## ----noisePvalues--------------------------------------------------------
noise_results <- pvaluesMultiEnrich(c("sample1", "sample2"), useGO, go_noise)

noise_pvalues <- pvalueDiffSig(noise_results$values, pCutoff=0.05, log=TRUE)

noise_pvalues$sizeClass <- useGOSizeClass
noise_pvalues$size <- goFracNoise[1:length(useGOSizeClass),4]
noise_pvalues$frac <- goFracNoise$frac[1:length(useGOSizeClass)]


## ----plotSummarizeNoise--------------------------------------------------
noise_pvalues$sizeClass <- factor(noise_pvalues$sizeClass, c("low", "med", "hi"), ordered=TRUE)
ggplot(noise_pvalues, aes(x=frac, y=diff, color=sigState)) + geom_point() + facet_grid(. ~ sizeClass, scales="free_x") + useTheme + ggtitle("Sample - Combined with Noise") + xlab("Fraction") + ylab("Diff. log p-values")


## ----summarizeCountsNoise------------------------------------------------
tapply(noise_pvalues$diff, noise_pvalues$sizeClass, function(x){sum(x > 0)})
tapply(noise_pvalues$diff, noise_pvalues$sizeClass, mean)


## ----diffDiffs-----------------------------------------------------------
diff_pvalues <- noise_pvalues
diff_pvalues$diff <- noise_pvalues$diff - noNoise_pvalues$diff
ggplot(diff_pvalues, aes(x=frac, y=diff)) + geom_point() + facet_grid(. ~ sizeClass, scales="free_x") + ggtitle("difference of differences noise vs noNoise") + useTheme + xlab("Fraction") + ylab("Diff. noise / noNoise")
tapply(diff_pvalues$diff, diff_pvalues$sizeClass, mean)


## ----geneSamples, eval=FALSE---------------------------------------------
## nTimes <- 100
## nGene <- 1000
## 
## testSamples <- lapply(seq(1, nTimes), function(x){
##   list(s1 = sampleTerms(useGO, hsGO, 1000, 4)[1:nGene],
##        s2 = sampleTerms(useGO, hsGO, 1000, 4)[1:nGene],
##        noise = sample(noiseGenes, nNoise))
## })
## 
## save(useGO, testSamples, universeGenes, useGOSizeClass, file="inst/data/100GeneSamples.RData")


## ----runSamples, eval=FALSE----------------------------------------------
## data("100GeneSamples")
## library(snowfall)
## # load("inst/data/100GeneSamples.RData")
## testFun <- function(x){
##   nnList <- x[c("s1", "s2")]
##   nnRes <- hyperGOMultiEnrichment(nnList, universe=universeGenes)
##   nn_res <- pvaluesMultiEnrich(c("s1", "s2"), useGO, nnRes)
##   nnP <- pvalueDiffSig(nn_res$values, pCutoff=0.05, log=TRUE)
##   nnP$class <- useGOSizeClass
## 
##   noiList <- list(s1 = c(x$s1, x$noise), s2 = c(x$s2, x$noise))
##   noiRes <- hyperGOMultiEnrichment(noiList, universe=universeGenes)
##   noi_res <- pvaluesMultiEnrich(c("s1", "s2"), useGO, noiRes)
##   noiP <- pvalueDiffSig(noi_res$values, pCutoff=0.05, log=TRUE)
##   noiP$class <- useGOSizeClass
## 
##   return(list(clean = nnP, noise = noiP))
## }
## sfInit(parallel=TRUE, cpus=10)
## sfExport("useGOSizeClass", "universeGenes", "useGO")
## sfLibrary(ccPaper)
## sfLibrary(GO.db)
## sfLibrary(org.Hs.eg.db)
## testRes <- sfLapply(testSamples, testFun)
## sfStop()
## save(testRes, file="inst/data/100GeneResults.RData")


## ----noNoiseStats--------------------------------------------------------
data("100GeneResults")
cleanRes <- lapply(testRes, function(x){x$clean})
testStats_clean <- sameGOStats(cleanRes)
testStats_clean$class <- useGOSizeClass
testStats_clean$frac <- goFractions$frac[1:length(useGOSizeClass)]
testStats_clean$class <- factor(testStats_clean$class, levels=c("low", "med", "hi"), ordered=TRUE)
p <- ggplot(testStats_clean, aes(x=frac, y=mean, ymax=max, ymin=min)) + geom_point() + geom_errorbar() + facet_grid(. ~ class, scales="free_x") + useTheme + ggtitle("Sample - List Differences, 100 Gene Samples") + xlab("Fraction") + ylab("Mean diff. log p-values")
p
savePlot(p, "genesamples100_sampleListDiffs.png")


## ----noiseStats----------------------------------------------------------
noiseRes <- lapply(testRes, function(x){x$noise})
testStats_noise <- sameGOStats(noiseRes)
testStats_noise$class <- useGOSizeClass
testStats_noise$frac <- goFracNoise$frac[1:length(useGOSizeClass)]
testStats_noise$class <- factor(testStats_noise$class, levels=c("low", "med", "hi"), ordered=TRUE)
ggplot(testStats_noise, aes(x=frac, y=mean, ymax=max, ymin=min)) + geom_point() + geom_errorbar() + facet_grid(. ~ class, scales="free_x") + ggtitle("Sample - List Differences, 100 Gene Samples, with noise") + xlab("Fraction") + ylab("Mean diff. log p-values") + useTheme


## ----100GeneNoiseClean---------------------------------------------------
noiseCleanDiff <- diffGOStats(testRes)
noiseCleanDiff$class <- useGOSizeClass
noiseCleanDiff$frac <- goFracNoise$frac[1:length(useGOSizeClass)]
noiseCleanDiff$class <- factor(noiseCleanDiff$class, levels=c("low", "med", "hi"), ordered=TRUE)
ggplot(noiseCleanDiff, aes(x=frac, y=mean, ymax=max, ymin=min)) + geom_point() + geom_errorbar() + facet_grid(. ~ class, scales="free_x") + ggtitle("Difference of Differences, 100 Gene Samples") + useTheme + xlab("Fraction") + ylab("Diff noise / noNoise")


## ----genSamples, eval=FALSE----------------------------------------------
## goLimits <- list(low = c(10, 100), med = c(250, 500), hi = c(500, 1500))
## 
## nSample <- 100
## 
## sfInit(parallel=TRUE, cpus=10)
## sfLibrary(ccPaper)
## sfExport("hsGO", "nGO", "goLimits")
## 
## testGOSample <- sfLapply(seq(1, nSample), function(x){
##   t1 <- goSample(hsGO, nGO, goLimits, nSample=2, nGene=1000, nNoise=500)
##   return(t1)
## })
## 
## sfStop()
## save(testGOSample, file="inst/data/100GOSamples.RData")


## ----runGOSamples, eval=FALSE--------------------------------------------
## data("100GOSamples")
## runGOSamples <- function(x){
##   cleanSample <- x$geneSample
##   cleanRes <- hyperGOMultiEnrichment(cleanSample, universe=universeGenes)
## 
##   cleanP <- pvaluesMultiEnrich(names(cleanSample), x$goData$id, cleanRes)
##   cleanP_dif <- pvalueDiffSig(cleanP$values, pCutoff=0.05, log=TRUE)
##   cleanP_dif <- cbind(cleanP_dif, x$goData)
## 
##   noiseSample <- lapply(names(cleanSample), function(inSample){
##     c(x$geneSample[[inSample]], x$noiseSample[[inSample]])
##   })
##   names(noiseSample) <- names(cleanSample)
##   noiseRes <- hyperGOMultiEnrichment(noiseSample, universe=universeGenes)
##   noiseP <- pvaluesMultiEnrich(names(cleanSample), x$goData$id, noiseRes)
##   noiseP_dif <- pvalueDiffSig(noiseP$values, pCutoff=0.05, log=TRUE)
##   noiseP_dif <- cbind(noiseP_dif, x$goData)
##   return(list(clean = cleanP_dif, noise = noiseP_dif))
## }
## 
## library(snowfall)
## sfInit(parallel=TRUE, cpus=10)
## sfLibrary(ccPaper)
## sfLibrary(org.Hs.eg.db)
## sfLibrary(GO.db)
## sfExport("universeGenes")
## 
## testGORes <- sfLapply(testGOSample, runGOSamples)
## sfStop()
## save(testGORes, file="inst/data/100GOResults.RData")


## ----lookDiffs-----------------------------------------------------------
data("100GOResults")

grab_noNoise <- function(x){
  x$clean[,c("diff", "class", "frac")]
}

go_noNoise <- lapply(testGORes, grab_noNoise)
go_noNoise <- do.call(rbind, go_noNoise)
go_noNoise$class <- factor(go_noNoise$class, levels=c("low", "med", "hi"), ordered=TRUE) 
go_noNoiseT <- go_noNoise[(!is.nan(go_noNoise$diff)),]
p <- ggplot(go_noNoiseT, aes(x = frac, y = diff)) + geom_point() + facet_grid(. ~ class, scales="free_x") + useTheme + stat_density2d(aes(fill = ..level..), geom="polygon") + ggtitle("Sample - List Differences, 100 GO Samples, noNoise") + xlab("Fraction") + ylab("Diff. log p-values") + labs(fill="density")
p
savePlot(p, "gosamples100_samplelistDiff.png")

grab_noise <- function(x){
  x$noise[, c("diff", "class", "frac")]
}
go_noise <- lapply(testGORes, grab_noise)
go_noise <- do.call(rbind, go_noise)
go_noise$class <- factor(go_noise$class, levels=c("low", "med", "hi"), ordered=TRUE)
ggplot(go_noise, aes(x = frac, y = diff)) + geom_point() + facet_grid(. ~ class, scales="free_x") + ggtitle("Sample - List Differences, 100 GO Samples, Noise") + xlab("Fraction") + ylab("Diff. log p-values") + useTheme

grab_noiseDiff <- function(x){
  outRes <- x$clean[, c("diff", "class", "frac")]
  outRes$diff <- x$noise$diff - x$clean$diff
  outRes
}
go_diff <- lapply(testGORes, grab_noiseDiff)
go_diff <- do.call(rbind, go_diff)
go_diff$class <- factor(go_diff$class, levels=c("low", "med", "hi"), ordered=TRUE)
ggplot(go_diff, aes(x = frac, y = diff)) + geom_point() + facet_grid(. ~ class, scales="free_x") + ggtitle("Difference of Differences, 100 GO Samples") + xlab("Fraction") + ylab("Diff. log p-values") + useTheme
ggplot(go_diff, aes(x = frac, y = diff, fill=class)) + geom_boxplot(position="identity", alpha=0.4) + ggtitle("Difference of Differences, 100 GO Samples") + xlab("Fraction") + ylab("Diff. noise / noNoise") + useTheme


## ----compNoiseNoNoise----------------------------------------------------
go_noise$noise <- "noise"
go_noNoise$noise <- "none"

go_all <- rbind(go_noise, go_noNoise)
rownames(go_all) <- NULL
ggplot(go_all, aes(x = frac, y = diff, fill = noise)) + geom_violin(trim=TRUE, alpha=0.5, position="identity") + facet_grid(. ~ class, scales="free_x") + ggtitle("Difference of Differences, 100 GO Samples") + xlab("Fraction") + ylab("Diff. log p-values") + useTheme


## ----aggregate-----------------------------------------------------------
go_diff$frac10 <- round(go_diff$frac * 10) / 10
ggplot(go_diff, aes(x = factor(frac10), y = diff)) + geom_boxplot() + facet_grid(. ~ class, scales="free_x") + ggtitle("Difference of Differences, 100 GO Samples") + xlab("Fraction") + ylab("Diff. log p-values") + useTheme

probPoint <- is.infinite(go_diff$diff) | is.na(go_diff$diff) | is.nan(go_diff$diff)
go_diff <- go_diff[!probPoint, ]
ggplot(go_diff, aes(x=frac, y=diff)) + stat_smooth() + facet_grid(. ~ class, scales="free_x") + ggtitle("Difference of Differences, 100 GO Samples") + xlab("Fraction") + ylab("Diff. log p-values") + useTheme


## ----sweepNoise, eval=FALSE----------------------------------------------
## noiseGenes <- possibleNoise(hsGO, useGO)
## sizeNoise <- seq(0, 1000, 10)
## fracShared <- seq(0, 1, 0.01)
## noiseSamples <- sweepNoiseSample(noiseGenes, nSamples=2, sizeNoise, fracShared)
## 
## nGene <- 1000
## cleanSample <- list(s1 = sampleTerms(useGO, hsGO, nGene, 4)[1:nGene],
##                     s2 = sampleTerms(useGO, hsGO, nGene, 4)[1:nGene])
## 
## save(noiseSamples, cleanSample, sizeNoise, fracShared, universeGenes, useGO, useGOSizeClass, file=file.path(useDir, 'inst/data/noiseSamples.RData'))


## ----runSweptNoise, eval=FALSE-------------------------------------------
## data(noiseSamples)
## 
## useNoise <- unlist(noiseSamples, recursive=FALSE)
## sizeNoiseAll <- rep(sizeNoise, each=length(fracShared))
## fracSharedAll <- rep(fracShared, length(sizeNoise))
## 
## runSweepNoise <- function(noiseGenes){
##   inNoise <- lapply(seq(1, length(cleanSample)), function(x){
##     c(cleanSample[[x]], noiseGenes[[x]])
##   })
##   names(inNoise) <- names(cleanSample)
## 
##   noiseRes <- hyperGOMultiEnrichment(inNoise, universe=universeGenes)
##   noiseP <- pvaluesMultiEnrich(names(inNoise), useGO, noiseRes)
##   noisePComp <- pvalueDiffSig(noiseP$values, pCutoff=0.05, log=TRUE)
##   noisePComp$class <- useGOSizeClass
##   return(noisePComp)
## }
## 
## library(snowfall)
## sfInit(parallel=TRUE, cpus=10)
## sfLibrary(ccPaper)
## sfLibrary(org.Hs.eg.db)
## sfLibrary(GO.db)
## sfExport("universeGenes", "useGO", "useGOSizeClass", "cleanSample")
## noiseSweepRes <- sfLapply(useNoise, runSweepNoise)
## sfStop()
## save(noiseSweepRes, file.path(useDir, "inst/data/noiseSweepRes.RData"))


## ----loadSweptNoise------------------------------------------------------
data(noiseSamples)
data(noiseSweepRes)

noiseFrac <- data.frame(size=rep(sizeNoise, each=length(fracShared)),
                        shared=rep(fracShared, length(sizeNoise)))


## ----examineSweptNoise---------------------------------------------------
outMedians <- lapply(noiseSweepRes, function(inRes){
  tapply(inRes$diff, inRes$class, median)
})
outMedians <- do.call(rbind, outMedians)
outMedians <- cbind(outMedians, noiseFrac)

ggplot(outMedians, aes(x=size, y=hi, color=shared)) + geom_point() + useTheme + ggtitle("Hi Median Difference") + xlab("Noise Added") + ylab("Median Diff.")

ggplot(outMedians, aes(x=size, y=med, color=shared)) + geom_point() + useTheme + ggtitle("Medium Median Difference") + xlab("Noise Added") + ylab("Median Diff.")

ggplot(outMedians, aes(x=size, y=low, color=shared)) + geom_point() + useTheme + ggtitle("Low Median Difference") + xlab("Noise Added") + ylab("Median Diff.")


## ----examine2------------------------------------------------------------
outMedians2 <- data.frame(diff=c(outMedians$low, outMedians$med, outMedians$hi), class=c(rep("low", nrow(outMedians)), rep("med", nrow(outMedians)), rep("hi", nrow(outMedians))), size=rep(outMedians$size, 3), shared=rep(outMedians$shared, 3))
outMedians2$class <- factor(outMedians2$class, c("low", "med", "hi"), ordered=TRUE)

p <- ggplot(outMedians2, aes(x=size, y=diff, color=shared)) + geom_point() + useTheme + facet_grid(class ~ ., scales="free_y") + ggtitle("Sample - Combined as a function of size and shared fraction")
p
savePlot(p, "samplecombined_functionSizeFraction.png")


## ----generateSampleRanks-------------------------------------------------
tALL <- topTable(fitALL2, number=Inf)
tALL <- tALL$t[(tALL$t < 10) & (tALL$t > -10)]
tValues <- runif(length(universeGenes), -6, 6)

gseaSamples <- list(s1 = sampleTerms(useGO, hsGO, 1000, 4),
                    s2 = sampleTerms(useGO, hsGO, 1000, 4))

gseaGOFracS1 <- calcFraction(hsGO[useGO], list(s1=gseaSamples[[1]]))

gseaTValues <- sampletValueGenGSEA(gseaSamples, universeGenes, tValues)


## ----checkValues---------------------------------------------------------
tmpDat <- data.frame(x=gseaTValues[[1]], type="all", stringsAsFactors=FALSE)
tmpDat$type[(universeGenes %in% gseaSamples[[1]])] <- "sample"
ggplot(tmpDat, aes(x=x, fill=type)) + geom_histogram(binwidth=0.1, position="identity") + ggtitle("location of samples") + useTheme

tmpComb <- data.frame(x=gseaTValues[["comb"]], type="all", stringsAsFactors=FALSE)
tmpComb$type[(universeGenes %in% gseaSamples[[1]])] <- "sample"
ggplot(tmpComb, aes(x=x, fill=type)) + geom_histogram(binwidth=0.1, position="identity") + ggtitle("location of samples") + useTheme


## ----unSomeGSEA----------------------------------------------------------
go2index <- symbols2indices(hsGO[useGO], universeGenes)
gseaSetResults <- lapply(gseaTValues, function(x){
  mmGeneSetTest(go2index, x)
})


## ----gseaBetterGO--------------------------------------------------------
cleanDiff <- gseaDiffs(gseaSetResults)
cleanDiff$sizeClass <- useGOSizeClass
cleanDiff$sizeClass <- factor(cleanDiff$sizeClass, c("low", "med", "hi"), ordered=TRUE)
cleanDiff$frac <- gseaGOFracS1$frac

ggplot(cleanDiff, aes(x=frac, y=diff)) + geom_point() + facet_grid(. ~ sizeClass) + ggtitle("Sample - List for GSEA, using GO") + useTheme + xlab("Fraction") + ylab("Diff. log p-values")


## ----gseaNoise-----------------------------------------------------------
nNoise <- 2000
noiseGenes <- possibleNoise(hsGO, useGO)
noiseGenes <- sample(noiseGenes, nNoise)

gseaNoise <- list(s1=c(gseaSamples[["s1"]], noiseGenes),
                  s2=c(gseaSamples[["s2"]], noiseGenes))

gseaTValues_noise <- sampletValueGenGSEA(gseaNoise, universeGenes, tValues)

tmpDat <- data.frame(x=gseaTValues_noise[[1]], type="all", stringsAsFactors=FALSE)
tmpDat$type[(universeGenes %in% gseaSamples[[1]])] <- "sample"
ggplot(tmpDat, aes(x=x, fill=type)) + geom_histogram(binwidth=0.1, position="identity", alpha=0.75) + ggtitle("Sample - List GSEA GO, adding noise genes") + useTheme


## ----runGSEANoise--------------------------------------------------------
gseaSetResults_noise <- lapply(gseaTValues_noise, function(x){
  mmGeneSetTest(go2index, x)
})


## ----gseaBetter----------------------------------------------------------
noiseDiff <- gseaDiffs(gseaSetResults_noise)
noiseDiff$sizeClass <- useGOSizeClass
noiseDiff$sizeClass <- factor(noiseDiff$sizeClass, c("low", "med", "hi"), ordered=TRUE)
noiseDiff$frac <- gseaGOFracS1$frac

ggplot(noiseDiff, aes(x=frac, y=diff)) + geom_point() + facet_grid(. ~ sizeClass) + ggtitle("Sample - List, GSEA GO with Noise") + useTheme + xlab("Fraction") + ylab("Diff. log p-values")

diffDiff <- data.frame(diff=(noiseDiff$diff - cleanDiff$diff), sizeClass=noiseDiff$sizeClass, frac=noiseDiff$frac)
ggplot(diffDiff, aes(x=frac, y=diff)) + geom_point() + facet_grid(. ~ sizeClass) + ggtitle("Difference of Differences, GSEA GO") + useTheme + xlab("Fraction") + ylab("Diff. noise / noNoise")


## ----independentGSEA-----------------------------------------------------
nGene <- 10000
indepTValues <- sort(runif(nGene, -6, 6))
names(indepTValues) <- as.character(seq(1, nGene))

gseaMap <- list(t1=sample(seq(1, nGene/2), 100))

tmpDat <- data.frame(tVal=indepTValues, status="all", stringsAsFactors=F)
tmpDat$status[gseaMap[["t1"]]] <- "sample"
ggplot(tmpDat, aes(x=tVal, fill=status)) + geom_histogram(binwidt=0.1, position="identity") + ggtitle("Check independent distribution, 100%") + useTheme + xlab("t-statistic")

geneSetTest(gseaMap$t1, indepTValues, alternative="down", ranks.only=FALSE)

gseaMap80 <- list(t1 = c(sample(seq(1, nGene/2), 80), 
                         sample(seq((nGene/2) + 1, nGene), 20 )))

geneSetTest(gseaMap80$t1, indepTValues, alternative="down", ranks.only=FALSE)

tmpDat <- data.frame(tVal=indepTValues, status="all", stringsAsFactors=F)
tmpDat$status[gseaMap80[["t1"]]] <- "sample"
ggplot(tmpDat, aes(x=tVal, fill=status)) + geom_histogram(binwidt=0.1, position="identity") + ggtitle("check independent distribution, 80%") + useTheme + xlab("t-statistic")


gseaMap50 <- list(t1 = c(sample(seq(1, nGene/2), 50),
                         sample(seq((nGene/2) + 1, nGene), 50)))

geneSetTest(gseaMap50$t1, indepTValues, alternative="down", ranks.only=FALSE)

gseaMapUnif <- list(t1 = sample(nGene, 100))
geneSetTest(gseaMapUnif$t1, indepTValues, alternative="down", ranks.only=FALSE)


## ----checkNabove---------------------------------------------------------
gseaAbove100 <- list(t1=sample(seq(100, nGene/2), 100))
geneSetTest(gseaAbove100$t1, indepTValues, alternative="down", ranks.only=FALSE)

gseaAbove500 <- list(t1=sample(seq(501, nGene/2), 100))
geneSetTest(gseaAbove500$t1, indepTValues, alternative="down", ranks.only=FALSE)

gseaAbove1000 <- list(t1=sample(seq(1001, nGene/2), 100))
geneSetTest(gseaAbove1000$t1, indepTValues, alternative="down", ranks.only=FALSE)

gseaAbove2000 <- list(t1=sample(seq(2001, nGene/2), 100))
geneSetTest(gseaAbove2000$t1, indepTValues, alternative="down", ranks.only=FALSE)

gseaAbove4000 <- list(t1=sample(seq(4001, nGene/2), 100))
geneSetTest(gseaAbove4000$t1, indepTValues, alternative="down", ranks.only=FALSE)


## ----gseaFunctions, eval=FALSE-------------------------------------------
## ?gseaSizeRange
## ?gseaProportion


## ----checkGSEAProp-------------------------------------------------------
testSize <- gseaSizeRange(c(500,500), 50)
useStatistics <- runif(10000, -6, 6)
useSig <- runif(50, .1, 1)
testStats <- gseaProportion(testSize, useSig, useStatistics)

fullIndex <- seq(1, length(useStatistics))
useCut <- 5000

realRanks <- lapply(testStats, function(inStat){
  
  apply(inStat$stats, 2, function(inCol){
    outRank <- rank(inCol)[inStat$index]
    sum(outRank < useCut) / length(inStat$index)
  })
  
})

realRanks <- do.call(rbind, realRanks)
plot(realRanks[,1], useSig)
plot(realRanks[,2], useSig)
plot(realRanks[,3], useSig)


## ----testFunctions-------------------------------------------------------
testSize <- gseaSizeRange(c(500,500), 50)
useStatistics <- runif(10000, -6, 6)
useSig <- runif(50, .1, 1)
testStats <- gseaProportion(testSize, useSig, useStatistics)

outP <- lapply(testStats, function(inStat){
  tmpP <- lapply(seq(1, ncol(inStat$stats)), function(useColumn){
    geneSetTest(inStat$index, inStat$stats[,useColumn], "down", ranks.only=FALSE)
  })
  tmpP <- do.call(cbind, tmpP)
  colnames(tmpP) <- colnames(inStat$stats)
  return(tmpP)
})

outP2 <- do.call(rbind, outP)
outP2 <- -1 * log(outP2)

tmpP <- data.frame(set=apply(outP2[, c("s1", "s2")], 1, min), list=outP2[,"comb"])
tmpP$diff <- tmpP$set - tmpP$list
tmpP$frac <- useSig
p <- ggplot(tmpP, aes(x=frac, y=diff)) + geom_point() + useTheme + ggtitle("Sample - Combined average t-statistics") + xlab("Fraction") + ylab("Diff. log p-values")
p
savePlot(p, "samplecombined_tstatistics.png")


## ----gseaVaryProp, eval=FALSE--------------------------------------------
## ?gseaProportionVary


## ----checkgseaProportionsVary--------------------------------------------
testSize <- gseaSizeRange(c(500,500), 50)
useStatistics <- runif(10000, -6, 6)
useSig <- list(s1=runif(50, .1, 1), s2=runif(50, .1, 1))
testStats <- gseaProportionVary(testSize, useSig, useStatistics)

fullIndex <- seq(1, length(useStatistics))
useCut <- 5000

realRanks <- lapply(testStats, function(inStat){
  
  apply(inStat$stats, 2, function(inCol){
    outRank <- rank(inCol)[inStat$index]
    sum(outRank < useCut) / length(inStat$index)
  })
  
})

realRanks <- do.call(rbind, realRanks)
plot(realRanks[,1], useSig[[1]])
plot(realRanks[,2], useSig[[2]])
plot(realRanks[,3], useSig[[1]])
plot(realRanks[,3], useSig[[2]])


## ----gseaProportions2Data------------------------------------------------
testSize <- gseaSizeRange(c(500,500), 50)
useStatistics <- runif(10000, -6, 6)
useSig <- list(s1=runif(50, .1, 1), s2=runif(50, .1, 1))
testStats <- gseaProportionVary(testSize, useSig, useStatistics)

outP <- lapply(testStats, function(inStat){
  tmpP <- lapply(seq(1, ncol(inStat$stats)), function(useColumn){
    geneSetTest(inStat$index, inStat$stats[,useColumn], "down", ranks.only=FALSE)
  })
  tmpP <- do.call(cbind, tmpP)
  colnames(tmpP) <- colnames(inStat$stats)
  return(tmpP)
})

outP2 <- do.call(rbind, outP)
outP2 <- -1 * log(outP2)

tmpP <- data.frame(set=apply(outP2[, c("s1", "s2")], 1, min), list=outP2[,"comb"])
tmpP$diff <- tmpP$set - tmpP$list

realRanks <- lapply(testStats, function(inStat){
  
  apply(inStat$stats, 2, function(inCol){
    outRank <- rank(inCol)[inStat$index]
    sum(outRank < useCut) / length(inStat$index)
  })
  
})

realRanks <- do.call(rbind, realRanks)

tmpP$frac <- realRanks[,3]
ggplot(tmpP, aes(x=frac, y=diff)) + geom_point() + ggtitle("Sample - Combined vary fraction") + useTheme + xlab("Fraction") + ylab("Diff. log p-values")


## ----generateRegression--------------------------------------------------
regData <- data.frame(y=c(-1, 0, 1), x=c(0, 0.5, 1))
regModel <- lm(y ~ x, regData)

testData <- runif(10000, 0, 1)
outData <- regModel$coefficients[1] + regModel$coefficients[2]*testData
hist(outData, 100)


## ----pvalRanks-----------------------------------------------------------
nGene <- 10000
indepTValues <- sort(runif(nGene, 0, 1))
names(indepTValues) <- as.character(seq(1, nGene))

gseaMap <- list(t1=sample(seq(1, nGene/4), 100))

tmpDat <- data.frame(tVal=indepTValues, status="all", stringsAsFactors=F)
tmpDat$status[gseaMap[["t1"]]] <- "sample"
ggplot(tmpDat, aes(x=tVal, fill=status)) + geom_histogram(binwidt=0.01, position="identity") + ggtitle("check p-values, 100%") + useTheme + xlab("P-values")

geneSetTest(gseaMap$t1, p2signed(regModel,indepTValues), alternative="down", ranks.only=FALSE)

gseaMap80 <- list(t1 = c(sample(seq(1, nGene/2), 80), 
                         sample(seq((nGene/2) + 1, nGene), 20 )))

geneSetTest(gseaMap80$t1, p2signed(regModel,indepTValues), alternative="down", ranks.only=FALSE)

tmpDat <- data.frame(tVal=indepTValues, status="all", stringsAsFactors=F)
tmpDat$status[gseaMap80[["t1"]]] <- "sample"
ggplot(tmpDat, aes(x=tVal, fill=status)) + geom_histogram(binwidt=0.01, position="identity") + ggtitle("check p-values, 80%") + useTheme + xlab("P-values")


gseaMap50 <- list(t1 = c(sample(seq(1, nGene/2), 50),
                         sample(seq((nGene/2) + 1, nGene), 50)))

geneSetTest(gseaMap50$t1, p2signed(regModel,indepTValues), alternative="down", ranks.only=FALSE)


## ----gseaPropFisher, eval=FALSE------------------------------------------
## ?gseaProportionFisher


## ----checkGSEAPropFisher-------------------------------------------------
testSize <- gseaSizeRange(c(500,500), 50)
useStatistics <- runif(10000, 0, 1)
useSig <- runif(50, .1, 1)
testStats <- gseaProportionFisher(testSize, useSig, useStatistics)

fullIndex <- seq(1, length(useStatistics))
useCut <- 5000

realRanks <- lapply(testStats, function(inStat){
  
  apply(inStat$stats, 2, function(inCol){
    outRank <- rank(inCol)[inStat$index]
    sum(outRank < useCut) / length(inStat$index)
  })
  
})

realRanks <- do.call(rbind, realRanks)
plot(realRanks[,1], useSig)
plot(realRanks[,2], useSig)
plot(realRanks[,3], useSig)


## ----testFishers---------------------------------------------------------
testSize <- gseaSizeRange(c(500,500), 50)
useStatistics <- runif(10000, 0, 1)
useSig <- runif(50, .1, 1)
testStats <- gseaProportionFisher(testSize, useSig, useStatistics)

outP <- lapply(testStats, function(inStat){
  tmpP <- lapply(seq(1, ncol(inStat$stats)), function(useColumn){
    geneSetTest(inStat$index, p2signed(regModel, inStat$stats[,useColumn]), "down", ranks.only=FALSE)
  })
  tmpP <- do.call(cbind, tmpP)
  colnames(tmpP) <- colnames(inStat$stats)
  return(tmpP)
})

outP2 <- do.call(rbind, outP)
outP2 <- -1 * log(outP2)

tmpP <- data.frame(set=apply(outP2[, c("s1", "s2")], 1, min), list=outP2[,"comb"])
tmpP$diff <- tmpP$set - tmpP$list
tmpP$frac <- useSig
ggplot(tmpP, aes(x=frac, y=diff)) + geom_point() + ggtitle("Sample - Combined Fishers") + useTheme + xlab("Fraction") + ylab("Diff. log p-values")

# try also using ranks only
outPRank <- lapply(testStats, function(inStat){
  tmpP <- lapply(seq(1, ncol(inStat$stats)), function(useColumn){
    geneSetTest(inStat$index, inStat$stats, ranks.only=TRUE)
  })
  tmpP <- do.call(cbind, tmpP)
  colnames(tmpP) <- colnames(inStat$stats)
  return(tmpP)
})

outP2Rank <- do.call(rbind, outPRank)
outP2Rank <- -1 * log(outP2Rank)

tmpP <- data.frame(set=apply(outP2Rank[, c("s1", "s2")], 1, min), list=outP2Rank[,"comb"])
tmpP$diff <- tmpP$set - tmpP$list
tmpP$frac <- useSig
ggplot(tmpP, aes(x=frac, y=diff)) + geom_point() + ggtitle("Sample - Combined Fishers Ranks") + useTheme + xlab("Fraction") + ylab("Diff. log p-values")


## ----gseaProportionMax---------------------------------------------------
?gseaProportionMax


## ----checkGSEAMax--------------------------------------------------------
testSize <- gseaSizeRange(c(500,500), 50)
useStatistics <- runif(10000, 0, 1)
useSig <- runif(50, .1, 1)
testStats <- gseaProportionMax(testSize, useSig, useStatistics)

fullIndex <- seq(1, length(useStatistics))
useCut <- 5000

realRanks <- lapply(testStats, function(inStat){
  
  apply(inStat$stats, 2, function(inCol){
    outRank <- rank(inCol)[inStat$index]
    sum(outRank < useCut) / length(inStat$index)
  })
  
})

realRanks <- do.call(rbind, realRanks)
plot(realRanks[,1], useSig)
plot(realRanks[,2], useSig)
plot(realRanks[,3], useSig)


## ----runGSEAMax----------------------------------------------------------
testSize <- gseaSizeRange(c(500,500), 50)
useStatistics <- runif(10000, 0, 1)
useSig <- runif(50, .1, 1)
testStats <- gseaProportionMax(testSize, useSig, useStatistics)

outP <- lapply(testStats, function(inStat){
  tmpP <- lapply(seq(1, ncol(inStat$stats)), function(useColumn){
    geneSetTest(inStat$index, p2signed(regModel, inStat$stats[,useColumn]), "down", ranks.only=FALSE)
  })
  tmpP <- do.call(cbind, tmpP)
  colnames(tmpP) <- colnames(inStat$stats)
  return(tmpP)
})

outP2 <- do.call(rbind, outP)
outP2 <- -1 * log(outP2)

tmpP <- data.frame(set=apply(outP2[, c("s1", "s2")], 1, min), list=outP2[,"comb"])
tmpP$diff <- tmpP$set - tmpP$list
tmpP$frac <- useSig
ggplot(tmpP, aes(x=frac, y=diff)) + geom_point() + ggtitle("Sample - Combined Max-P") + useTheme + xlab("Fraction") + ylab("Diff. log p-values")


## ----varySizesAndProps---------------------------------------------------
gseaLo <- gseaSizeRange(c(20, 250), 50)
gseaMed <- gseaSizeRange(c(250, 500), 30)
gseaHi <- gseaSizeRange(c(500, 1000), 20)

gseaSizes <- c(gseaLo, gseaMed, gseaHi)
gseaSizeClass <- c(rep("low", 50), rep("med", 30), rep("hi", 20))

gseaFrac <- runif(100, 0.3, 1)
gseaStatistics <- runif(10000, 0, 1)

testStats <- gseaProportionMax(gseaSizes, gseaFrac, gseaStatistics)

outP <- lapply(testStats, function(inStat){
  tmpP <- lapply(seq(1, ncol(inStat$stats)), function(useColumn){
    geneSetTest(inStat$index, p2signed(regModel, inStat$stats[,useColumn]), "down", ranks.only=FALSE)
  })
  tmpP <- do.call(cbind, tmpP)
  colnames(tmpP) <- colnames(inStat$stats)
  return(tmpP)
})

outP2 <- do.call(rbind, outP)
outP2 <- -1 * log(outP2)

tmpP <- data.frame(set=apply(outP2[, c("s1", "s2")], 1, min), list=outP2[,"comb"])
tmpP$diff <- tmpP$set - tmpP$list
tmpP$frac <- gseaFrac
tmpP$class <- gseaSizeClass
tmpP$class <- factor(tmpP$class, c("low", "med", "hi"), ordered=TRUE)
ggplot(tmpP, aes(x=frac, y=diff)) + geom_point() + facet_grid(. ~ class) + useTheme

p <- ggplot(tmpP, aes(x=frac, y=diff)) + geom_point() + useTheme + ggtitle("Sample - Combined Max-P") + xlab("Fraction") + ylab("Diff. log p-values")
p
savePlot(p, "samplecombined_maxp.png")


## ----maxPRankOnly--------------------------------------------------------
outPRank <- lapply(testStats, function(inStat){
  tmpP <- lapply(seq(1, ncol(inStat$stats)), function(useColumn){
    geneSetTest(inStat$index, inStat$stats[,useColumn], ranks.only=TRUE)
  })
  tmpP <- do.call(cbind, tmpP)
  colnames(tmpP) <- colnames(inStat$stats)
  return(tmpP)
})

outP2Rank <- do.call(rbind, outPRank)
outP2Rank <- -1 * log(outP2Rank)

tmpPRank <- data.frame(set=apply(outP2Rank[, c("s1", "s2")], 1, min), list=outP2Rank[,"comb"])
tmpPRank$diff <- tmpPRank$set - tmpPRank$list
tmpPRank$frac <- gseaFrac
tmpPRank$class <- gseaSizeClass
tmpPRank$class <- factor(tmpPRank$class, c("low", "med", "hi"), ordered=TRUE)
ggplot(tmpPRank, aes(x=frac, y=diff)) + geom_point() + facet_grid(. ~ class) + ggtitle("Sample - Combined Max-P Ranks") + xlab("Fraction") + ylab("Diff. log p-values")


