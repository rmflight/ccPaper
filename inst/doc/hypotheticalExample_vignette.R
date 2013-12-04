
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


## ----checkFractions------------------------------------------------------
goFractions <- calcFraction(hsGO[useGO], list(sample1=sample1_org, sample2=sample2_org))
ggplot(goFractions, aes(x=size, y=frac, color=genelist)) + geom_point()


## ----runCalcs------------------------------------------------------------
samples_noNoise <- list(sample1 = sample1_org, sample2 = sample2_org)
go_noNoise <- hyperGOMultiEnrichment(samples_noNoise, universe=universeGenes)


## ----noNoisePvalues------------------------------------------------------
noNoise_pvalues <- pvaluesMultiEnrich(c("sample1", "sample2"), useGO, go_noNoise)

noNoise_pvalues <- pvalueDiffSig(noNoise_pvalues, pCutoff=0.05, log=TRUE)


## ----noNoiseAddInfo------------------------------------------------------
sizeClass <- rep(c("low", "med", "hi"), each=20)
names(sizeClass) <- useGO
noNoise_pvalues$sizeClass <- sizeClass
noNoise_pvalues$size <- goFractions[1:60,4]
noNoise_pvalues$frac <- goFractions$frac[1:60]


## ----noNoisePlotStuff----------------------------------------------------
noNoise_pvalues$sizeClass <- factor(noNoise_pvalues$sizeClass, levels=c("low", "med", "hi"), ordered=TRUE)
ggplot(noNoise_pvalues, aes(x=frac, y=diff, color=sigState)) + geom_point() + facet_grid(. ~ sizeClass, scales="free_x")
ggplot(noNoise_pvalues, aes(x=sizeClass, y=diff)) + geom_boxplot() + geom_point()


## ----noNoiseBetterPVals--------------------------------------------------
tapply(noNoise_pvalues$diff, noNoise_pvalues$sizeClass, function(x){sum(x > 0)})


## ----noNoiseBetterMean---------------------------------------------------
tapply(noNoise_pvalues$diff, noNoise_pvalues$sizeClass, mean)


## ----noiseGenes----------------------------------------------------------
nNoise <- 500
gene2HsGO <- reverseSplit(hsGO)
not_useGO <- sapply(gene2HsGO, function(x){
  sum(x %in% useGO) == 0
})
noiseGenes <- names(not_useGO)[not_useGO]
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
noise_pvalues <- pvaluesMultiEnrich(c("sample1", "sample2"), useGO, go_noise)

noise_pvalues <- pvalueDiffSig(noise_pvalues, pCutoff=0.05, log=TRUE)

noise_pvalues$sizeClass <- sizeClass
noise_pvalues$size <- goFracNoise[1:60,4]
noise_pvalues$frac <- goFracNoise$frac[1:60]


## ----plotSummarizeNoise--------------------------------------------------
noise_pvalues$sizeClass <- factor(noise_pvalues$sizeClass, c("low", "med", "hi"), ordered=TRUE)
ggplot(noise_pvalues, aes(x=frac, y=diff, color=sigState)) + geom_point() + facet_grid(. ~ sizeClass, scales="free_x")


## ----summarizeCountsNoise------------------------------------------------
tapply(noise_pvalues$diff, noise_pvalues$sizeClass, function(x){sum(x > 0)})
tapply(noise_pvalues$diff, noise_pvalues$sizeClass, mean)


