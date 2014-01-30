
## ----customCSS, include=FALSE--------------------------------------------
cssFile <- system.file("extdata", "style.css", package="ccPaperRev")
options(markdown.HTML.stylesheet = cssFile)


## ----loadLibrary, message=FALSE------------------------------------------
library(ccPaperRev)


## ----loadLibraries-------------------------------------------------------
options(stringsAsFactors=TRUE)
library(affy)


## ----skinRaw, eval=FALSE-------------------------------------------------
## currLoc <- getwd()
## dataLoc <- "/mlab/data/rmflight/Documents/projects/work/Petruska/drg_sc_skin_denervation/Skin/"
## setwd(dataLoc)
## 
## attFile <- "Sample information file.txt"
## dataAttributes <- read.table(attFile, sep="\t", header=T)
## dataAttributes$timePoint <- 0
## dataAttributes$timePoint[(dataAttributes$Time == "7 day")] <- 7
## dataAttributes$timePoint[(dataAttributes$Time == "14 day")] <- 14
## 
## celData <- ReadAffy(phenoData=dataAttributes)
## skinData <- rma(celData)
## 
## setwd(currLoc)
## .sessionInfo <- sessionInfo()
## .timeDate <- Sys.time()
## save(skinData, .sessionInfo, .timeDate, file="inst/data/skin_rma_data.RData")


## ----muscleRaw, eval=FALSE-----------------------------------------------
## library(GEOquery)
## muscleData <- getGEO(GEO="GSE4411")[[1]]
## 
## controlSamples <- grep("Control", pData(muscleData)$title)
## muscleData <- muscleData[,controlSamples]
## 
## tmpP <- pData(muscleData)
## tmpP$innervation <- "innervated"
## denLoc <- grep("Denervated", tmpP$title)
## tmpP$innervation[denLoc] <- "denervated"
## pData(muscleData) <- tmpP
## .sessionInfo <- sessionInfo()
## .timeDate <- Sys.time()
## save(muscleData, .sessionInfo, .timeDate, file="inst/data/muscle_data.RData")


## ----diffSkin------------------------------------------------------------
library(limma)
library(rat2302.db)
data(skin_rma_data)

skinComps <- c("T7 - T0", "T14 - T0")
geneID <- unlist(mget(featureNames(skinData), rat2302ENTREZID))

skinExpr <- exprs(skinData)
skinCollapse <- collapseProbes(skinExpr, geneID) # collapse to single genes using median of expression
skinCharacter <- pData(skinData)
skinCharacter$timePoint2 <- paste("T", skinCharacter$timePoint, sep="")

skinFC <- rankGenes(skinCollapse, skinCharacter$timePoint2, skinComps, doAggregation=FALSE, aggregateIndex=NA)
names(skinFC) <- c("T7", "T14")

skinDiff <- lapply(skinFC, getDiffGenes, id="id")


## ----setupGeneSets-------------------------------------------------------
library(GO.db)
library(org.Rn.eg.db)
library(limma)
rnGO <- as.list(org.Rn.egGO2ALLEGS)
rnGO <- rnGO[(Ontology(names(rnGO)) == "BP")]
rnGO <- lapply(rnGO, unique)

skinSet <- symbols2indices(rnGO, rownames(skinCollapse))


## ----skinContrasts, eval=FALSE-------------------------------------------
## skinStatus <- skinCharacter$timePoint2
## f <- factor(skinStatus)
## skinDesign <- model.matrix(~0 + f)
## colnames(skinDesign) <- levels(f)
## 
## skinContrast <- makeContrasts(contrasts=skinComps, levels=skinDesign)
## 
## save(list=ls(), file="runSkinRomer.RData")


## ----skinRomer, eval=FALSE-----------------------------------------------
## load("runSkinRomer.RData")
## options(mc.cores=12) # at least if we are on hera
## t1 <- Sys.time()
## skinRomer <- multicontrastRomer(skinSet, skinCollapse, skinDesign, skinContrast, nrot=10000)
## names(skinRomer) <- c("T7", "T14")
## save(skinRomer, file="inst/data/skinRomer.RData")
## t2 <- Sys.time()
## difftime(t2, t1)


## ----runGeneSetTest, eval=FALSE------------------------------------------
## load("runSkinRomer.RData")
## options(mc.cores=12) # at least if we are on hera
## t1 <- Sys.time()
## skinGeneSetTest <- lapply(skinFC, function(x){
##   mmGeneSetTest(skinSet, x$t)
## })
## 
## t2 <- Sys.time()
## difftime(t2, t1)
## save(skinGeneSetTest, t1, t2, file="inst/data/skinGeneSetTest.RData")


## ----diffMuscle----------------------------------------------------------
library(limma)
library(mouse4302.db)
data(muscle_data)

muscleComps <- c("denervated - innervated")
geneID <- unlist(mget(featureNames(muscleData), mouse4302ENTREZID))

muscleExpr <- exprs(muscleData)
muscleCollapse <- collapseProbes(muscleExpr, geneID) # collapse to single genes using median of expression
muscleCharacter <- pData(muscleData)

muscleFC <- rankGenes(muscleCollapse, muscleCharacter$innervation, muscleComps, doAggregation=FALSE, aggregateIndex=NA)
names(muscleFC) <- "denervation"

muscleDiff <- lapply(muscleFC, getDiffGenes, id="id")


## ----setupGeneSetsMuscle-------------------------------------------------
library(GO.db)
library(org.Mm.eg.db)
library(limma)
mmGO <- as.list(org.Mm.egGO2ALLEGS)
mmGO <- mmGO[(Ontology(names(mmGO)) == "BP")]
mmGO <- lapply(mmGO, unique)

muscleSet <- symbols2indices(mmGO, rownames(muscleCollapse))


## ----muscleContrasts, eval=FALSE-----------------------------------------
## muscleStatus <- muscleCharacter$innervation
## f <- factor(muscleStatus)
## muscleDesign <- model.matrix(~0 + f)
## colnames(muscleDesign) <- levels(f)
## 
## muscleContrast <- makeContrasts(contrasts=muscleComps, levels=muscleDesign)
## 
## save(muscleSet, muscleCollapse, muscleDesign, muscleContrast, mmGO, muscleFC, file="runMuscleRomer.RData")


## ----muscleRomer, eval=FALSE---------------------------------------------
## load("runMuscleRomer.RData")
## options(mc.cores=12) # at least if we are on hera
## t1 <- Sys.time()
## muscleRomer <- multicontrastRomer(muscleSet, muscleCollapse, muscleDesign, muscleContrast, nrot=10000)
## names(muscleRomer) <- "denervation"
## save(muscleRomer, file="inst/data/muscleRomer.RData")
## t2 <- Sys.time()
## difftime(t2, t1)


## ----getDataOut----------------------------------------------------------
data(muscleRomer)
data(skinRomer)

bothGO <- union(names(rnGO), names(mmGO))
combMapping <- lapply(bothGO, function(inName){
  isMM <- inName %in% names(mmGO)
  isRN <- inName %in% names(rnGO)
  
  tmpGene <- character(0)
  
  if (isMM){
    tmpGene <- c(tmpGene, mmGO[[inName]])
  }
  if (isRN){
    tmpGene <- c(tmpGene, rnGO[[inName]])
  }
  tmpGene <- unique(tmpGene)
  return(tmpGene)
})

names(combMapping) <- bothGO

geneAnnMapping <- new("namedList", .Data=combMapping, names=names(combMapping))


getSigID <- function(inRomer, pCut=0.05, whichCol=c("Up", "Down")){
  sigID <- lapply(whichCol, function(inCol){
    rownames(inRomer[(inRomer[, inCol] <= pCut),])
  })
  names(sigID) <- whichCol
  return(sigID)
}

muscleSigRomer <- lapply(muscleRomer, getSigID, pCut=0.01)
muscleSigRomer <- unlist(muscleSigRomer, recursive=FALSE)
names(muscleSigRomer) <- c("muscle.Up", "muscle.Down")

skinSigRomer <- lapply(skinRomer, getSigID, pCut = 0.01)
skinSigRomer <- unlist(skinSigRomer, recursive=FALSE)

allSigRomer <- c(skinSigRomer, muscleSigRomer)

genCCSigList <- function(inSig){
  tmp <- new("ccSigList", sigID=inSig)
}

allCCSig <- lapply(allSigRomer, genCCSigList)

allCCRomer <- new("GENccEnrichResult", allCCSig, categoryName="GSEAGO", geneAnnMapping=geneAnnMapping, overlapType="overlap", annDescription=Term(names(geneAnnMapping)))

allRomerOpts <- new("ccOptions", listNames=names(allCCRomer), colorType="pie")
compareSMRomer <- ccCompare(allCCRomer, allRomerOpts)
compareSMRomer


## ----saveData, eval=FALSE------------------------------------------------
## save(compareSMRomer, allRomerOpts, file="skinMuscleCompare.RData")


## ----out2cytoscape, eval=FALSE-------------------------------------------
## compareSMRomer <- breakEdges(compareSMRomer, 0.8)
## cwSM <- ccOutCyt(compareSMRomer, allRomerOpts, postText="GSEAGO", rpcPort=9001)


## ----runGeneSetTestMuscle, eval=FALSE------------------------------------
## load("runMuscleRomer.RData")
## options(mc.cores=12) # at least if we are on hera
## t1 <- Sys.time()
## muscleGeneSetTest <- lapply(muscleFC, function(x){
##   mmGeneSetTest(muscleSet, x$t)
## })
## 
## t2 <- Sys.time()
## difftime(t2, t1)
## save(muscleGeneSetTest, t1, t2, file="inst/data/muscleGeneSetTest.RData")


