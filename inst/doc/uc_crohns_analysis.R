
## ----customCSS, include=FALSE--------------------------------------------
cssFile <- system.file("extdata", "style.css", package="ccPaper")
options(markdown.HTML.stylesheet = cssFile)


## ----loadLibrary, message=FALSE------------------------------------------
library(ccPaper)


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


## ----saveNodeGroups, eval=FALSE------------------------------------------
## cucNodeGroups <- cytOutNodes("response to lipopolysaccharide and bacterial - UC.Up", cucRomerCW)
## cucNodeGroups <- cytOutNodes("regulation of inflammatory response - UC.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("hydrogen peroxide metabolism - UC.Up,CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("regulation of cell cycle and DNA damage response - UC.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("regulation of ubiquitination and ligase activity - UC.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("nucleotide and nucleoside metabolism - UC.Up,CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("amine metabolism - UC.Up,CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("glandular cell differentiation - UC.Down", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("oligodendrocyte differentiation - UC.Down,CROHNS.Down", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("cellular pattern specification - UC.Down,CROHNS.Down", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("NAD biosynthesis - CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("hormone metabolism - CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("extrinsic signal transduction - UC.Up,CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("regulation of nitric-oxide synthase - UC.Up,CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("fatty-acyl-CoA biosynthesis - UC.Up,CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("response to growth hormone - CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("melanin metabolism - CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("protein dephosphorylation - CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("chemokine and cytokine production - UC.Up,CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("membrane biogenesis and assembly - UC.Down", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("activin receptor signaling - UC.Down", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("nik/nk-kappab cascade - UC.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("er unfolded protein response - UC.Up, CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("regulation of ras/rac/rho gtpase activity - UC.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("antigen processing and presentation - UC.UP,CROHNS.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("COPII vesicle coating and targeting - UC.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("negative regulation of peptidase activity - UC.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("response to type 1 interferon - UC.Up", cucRomerCW, cucNodeGroups)
## cucNodeGroups <- cytOutNodes("protein N-linked glycosylation - UC.Up", cucRomerCW, cucNodeGroups)
## .sessionInfo <- sessionInfo()
## .timeDate <- Sys.time()
## save(cucRomerCW, cucNodeGroups, compareCUCRomer, cucRomerOpts, .sessionInfo, .timeDate, file="inst/data/cucCCOutput.RData")


## ----writeTableResults---------------------------------------------------
data(cucCCOutput)
allDescStrings <- sapply(cucNodeGroups, function(x){x$descStr})
string2List <- strsplit(allDescStrings, " - ", fixed=TRUE)
justDesc <- sapply(string2List, function(x){x[1]})
listMem <- sapply(string2List, function(x){x[2]})
listMemSplit <- strsplit(listMem, ",", fixed=TRUE)

descMembershipTable <- matrix("", nrow=length(allDescStrings), ncol=5)
colnames(descMembershipTable) <- c("Description", "UC.Down", "UC.Up", "CROHNS.Down", "CROHNS.Up")

descMembershipTable[,"Description"] <- justDesc

for (inRow in seq(1, nrow(descMembershipTable))){
  useSplit <- listMemSplit[[inRow]]
  trimSplit <- gsub(" ", "", useSplit)
  useLocs <- sapply(trimSplit, grep, colnames(descMembershipTable), ignore.case=TRUE)
  descMembershipTable[inRow, useLocs] <- "X"
}

orderBy <- c("UC.Up,CROHNS.Up",
             "UC.UP,CROHNS.Up",
             "UC.Up",
             "UC.Down",
             "UC.Down,CROHNS.Down",
             "CROHNS.Up",
             "CROHNS.Down")

listMem <- gsub(" ", "", listMem)
newOrder <- unlist(lapply(orderBy, function(x){
  which(listMem %in% x)
}))

descMembershipTable <- descMembershipTable[newOrder,]
require(xtable)

# add an html link to each entry in the table
useLink <- paste('<a href="#loc', seq(1, nrow(descMembershipTable)), '">', descMembershipTable[,"Description"], '</a>', sep="")
descMembershipTable[, "Description"] <- useLink


## ----printTable, echo=FALSE, results='asis'------------------------------
cat('<a name="tableLink"></a>\n') # a link back to the table
print(xtable(descMembershipTable), type = 'html', html.table.attributes = 'style="border-spacing:20px 5px;"', include.rownames=FALSE, sanitize.text.function=function(x){x})


## ----detailTable, echo=FALSE, results='asis'-----------------------------
cucNodeGroups <- cucNodeGroups[newOrder]

for (iNode in seq(1:length(cucNodeGroups))){
  nodeSet <- cucNodeGroups[[iNode]]
  toLink <- paste('<a name="loc', iNode, '"></a>', sep="")
  cat(toLink, "\n")
  cat("####", nodeSet$descStr, sep=" ")
  cat("\n")
  
  nGO <- length(nodeSet$nodes)
  outMatrix <- matrix("", nrow=nGO, ncol=6)
  colnames(outMatrix) <- c("GOID", "Description", "UC.Down", "UC.Up", "CROHNS.Down", "CROHNS.Up")
  rownames(outMatrix) <- nodeSet$nodes
  outMatrix[,"GOID"] <- nodeSet$nodes
  
  for (useGO in nodeSet$nodes){
    outMatrix[useGO, "Description"] <- nodeSet$nodeData[[useGO]]$Desc
    listMems <- strsplit(nodeSet$nodeData[[useGO]]$listMembership, ",")[[1]]
    listMems <- listMems[(nchar(listMems) > 0)]
    outMatrix[useGO, listMems] <- "X"
  }
  allMems <- sapply(nodeSet$nodeData, function(x){x$listMembership})
  reOrder <- order(allMems, decreasing=TRUE)
  outMatrix <- outMatrix[reOrder,]
  print(xtable(outMatrix), type = 'html', html.table.attributes = 'style="border-spacing:20px 5px;"', include.rownames=FALSE)
  cat('\n<a href="#tableLink">back</a>\n')
}


## ----getAllData----------------------------------------------------------
library(graph)
data(cucRomer)
data(cucCCOutput)

getTable <- function(inName){
  tmpDat <- cucRomer[[inName]]
  colnames(tmpDat) <- paste(inName, colnames(tmpDat), sep=".")
  tmpDat
}

cucTables <- lapply(names(cucRomer), getTable)
all.equal(rownames(cucTables[[1]]), rownames(cucTables[[2]]))

cucTables <- do.call(cbind, cucTables)
cucTables <- as.data.frame(cucTables, stringsAsFactors=FALSE)
cucTables$ID <- rownames(cucTables)

sigNodes <- nodes(compareCUCRomer@mainGraph)
nodeMembership <- unlist(nodeData(compareCUCRomer@mainGraph, sigNodes, "listMembership"))

cucTables <- cucTables[sigNodes,]
cucTables$membership <- nodeMembership
cucTables$description <- unlist(nodeData(compareCUCRomer@mainGraph, sigNodes, "Desc"))

groupLabel <- "nucleotide and nucleoside metabolism - UC.Up,CROHNS.Up"
allLabel <- sapply(cucNodeGroups, function(x){x$descStr})
whichLabel <- which(allLabel %in% groupLabel)
groupNodes <- cucNodeGroups[[whichLabel]]$nodes
cucTables[groupNodes,]


## ----checkRanks----------------------------------------------------------
newOrder <- order(cucTables[, 'UC.Up'], cucTables[, 'CROHNS.Up'], cucTables[, 'membership'])
cucTables <- cucTables[newOrder,]
which(cucTables$ID %in% groupNodes)


## ------------------------------------------------------------------------
Sys.time()
sessionInfo()


