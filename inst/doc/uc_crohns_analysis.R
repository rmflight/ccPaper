
## ----customCSS, include=FALSE--------------------------------------------
cssFile <- system.file("extdata", "style.css", package="ccPaperRev")
options(markdown.HTML.stylesheet = cssFile)


## ----loadLibrary, message=FALSE------------------------------------------
library(ccPaperRev)


## ----getData, eval=FALSE-------------------------------------------------
## library(GEOquery)
## 
## cucData <- getGEO(GEO="GSE36807")[[1]]
## 
## tmpP <- pData(cucData)
## tmpP$status <- "control"
## 
## cLoc <- grep("Crohn", tmpP$characteristics_ch1)
## tmpP$status[cLoc] <- "crohns"
## 
## ucLoc <- grep("Ulcerative colitis", tmpP$characteristics_ch1)
## tmpP$status[ucLoc] <- "uc"
## pData(cucData) <- tmpP
## 
## .sessionInfo <- sessionInfo()
## .timeDate <- Sys.time()
## save(cucData, .sessionInfo, .timeDate, file="inst/data/uc_crohns_rawData.RData")


## ----rankProbes----------------------------------------------------------
library(limma)
library(hgu133plus2.db)
data(uc_crohns_rawData)

cucComps <- c("uc - control", "crohns - control")
geneID <- unlist(mget(featureNames(cucData), hgu133plus2ENTREZID))

cucExpr <- exprs(cucData)
cucCollapse <- collapseProbes(cucExpr, geneID) # collapse to single genes using median of expression
cucCharacter <- pData(cucData)

cucFC <- rankGenes(cucCollapse, cucCharacter$status, cucComps, doAggregation=FALSE, aggregateIndex=NA)
names(cucFC) <- c("UC", "CROHNS")

cucDiff <- lapply(cucFC, getDiffGenes, id="id")
# this object can be run using base ccEnrich, or can be easily coerced


## ----ORA-----------------------------------------------------------------
library(GO.db)
cucGeneList <- list(UC_up=list(genes=cucDiff$UC$up, universe=cucDiff$UC$universe, annotation="org.Hs.eg.db"),
                    UC_dn=list(genes=cucDiff$UC$dn, universe=cucDiff$UC$universe, annotation="org.Hs.eg.db"),
                    CR_up=list(genes=cucDiff$CROHNS$up, universe=cucDiff$CROHNS$universe, annotation="org.Hs.eg.db"),
                    CR_dn=list(genes=cucDiff$CROHNS$dn, universe=cucDiff$CROHNS$universe, annotation="org.Hs.eg.db"))
cucGeneList <- new("ccGeneList", cucGeneList, ccType="BP")

cucEnrich <- ccEnrich(cucGeneList)
pvalueType(cucEnrich) <- "pval"
pvalueCutoff(cucEnrich) <- 0.01

cucEnrich

cucOpts <- new("ccOptions", listNames=names(cucGeneList), compareNames=
                c("UC_up", "UC_dn", "CR_up", "CR_dn", "UC_up,CR_up", "UC_up,CR_dn", "UC_dn,CR_up", "UC_dn,CR_dn"))
cucCompare <- ccCompare(cucEnrich, cucOpts)

cucCompare



## ----out2Cytoscape, eval=FALSE-------------------------------------------
## cucCompare <- breakEdges(cucCompare$BP, 0.8)
## cucCy <- ccOutCyt(cucCompare, cucOpts)
## breakEdges(cucCy, 1)


## ----setupGeneSets-------------------------------------------------------
library(GO.db)
library(org.Hs.eg.db)
library(limma)
hsGO <- as.list(org.Hs.egGO2ALLEGS)
hsGO <- hsGO[(Ontology(names(hsGO)) == "BP")]
hsGO <- lapply(hsGO, unique)

useSet <- symbols2indices(hsGO, rownames(cucCollapse))


## ----setupContrasts------------------------------------------------------
sampleStatus <- cucCharacter$status
doComps <- cucComps
f <- factor(sampleStatus)
design <- model.matrix(~0 + f)
colnames(design) <- levels(f)

contrast.matrix <- makeContrasts(contrasts=doComps, levels=design)


## ----runRomer, eval=FALSE------------------------------------------------
## cucRomer <- multicontrastRomer(useSet, cucCollapse, design, contrast.matrix, nrot=10000)
## names(cucRomer) <- c("UC", "CROHNS")
## save(cucRomer, file="inst/data/cucRomer.RData")
## 


## ----examineResults------------------------------------------------------
data(cucRomer)

geneAnnMapping <- new("namedList", .Data=hsGO, names=names(hsGO))

getSigID <- function(inRomer, pCut=0.05, whichCol=c("Up", "Down")){
  sigID <- lapply(whichCol, function(inCol){
    rownames(inRomer[(inRomer[, inCol] <= pCut),])
  })
  names(sigID) <- whichCol
  return(sigID)
}

cucSigRomer <- lapply(cucRomer, getSigID, pCut=0.01)
cucSigRomer <- unlist(cucSigRomer, recursive=FALSE)

genCCSigList <- function(inSig){
  tmp <- new("ccSigList", sigID=inSig)
}

cucCCSig <- lapply(cucSigRomer, genCCSigList)

cucCCRomer <- new("GENccEnrichResult", cucCCSig, categoryName="GSEAGO", geneAnnMapping=geneAnnMapping, overlapType="overlap", annDescription=Term(names(geneAnnMapping)))

cucRomerOpts <- new("ccOptions", listNames=names(cucCCRomer), compareNames=c(
  "UC.Up", "UC.Down", "CROHNS.Up", "CROHNS.Down", "UC.Up,CROHNS.Up", "UC.Up,CROHNS.Down", "UC.Down,CROHNS.Up", "UC.Down,CROHNS.Down"))
compareCUCRomer <- ccCompare(cucCCRomer, cucRomerOpts)
compareCUCRomer


## ----send2Cytoscape, eval=FALSE------------------------------------------
## compareCUCRomer <- breakEdges(compareCUCRomer, 0.8)
## cucRomerCW <- ccOutCyt(compareCUCRomer, cucRomerOpts, postText="romer", rpcPort=9001)
## generateLegend(cucRomerOpts)


