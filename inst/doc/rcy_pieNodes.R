
## ----customCSS, include=FALSE--------------------------------------------
cssFile <- system.file("extdata", "style.css", package="ccPaperRev")
options(markdown.HTML.stylesheet = cssFile)


## ----exampleGraph, eval=FALSE--------------------------------------------
## library(RCytoscape)
## g <- new ('graphNEL', edgemode='directed')
## g <- graph::addNode ('A', g)
## g <- graph::addNode ('B', g)
## g <- graph::addNode ('C', g)
## cw <- new.CytoscapeWindow ('vignette', graph=g)
## displayGraph(cw)
## layoutNetwork(cw, layout.name="grid")
## redraw(cw)


## ----pieGraphs-----------------------------------------------------------
pieData <- list(
  A = list(
    area = c(0.5, 0.5),
    color = c("blue", "green")),
  B = list(
    area = c(0.2, 0.8),
    color = c("red", "yellow")),
  C = list(
    area = c(0.6, 0.4),
    color = c("green", "purple"))
  )

plotPie <- function(inData){
  names(inData$area) <- ""
  pie(inData$area, col=inData$color, clockwise=TRUE)
}

plotPie(pieData[["A"]])
plotPie(pieData[["B"]])
plotPie(pieData[["C"]])


## ----savePie, eval=FALSE-------------------------------------------------
## saveDir <- file.path(getwd(), "vignettes/")
## 
## outFiles <- sapply(names(pieData), function(pieName){
##   pieFile <- paste("pieFile_", pieName, ".png", sep="")
##   outFile <- file.path(saveDir, pieFile)
##   png(outFile, bg="transparent")
##   plotPie(pieData[[pieName]])
##   dev.off()
##   return(outFile)
## })
## names(outFiles) <- NULL
## outFiles2 <- paste("file://localhost/", outFiles, sep="")


## ----showInCy, eval=FALSE------------------------------------------------
## setNodeImageDirect(cw, names(pieData), outFiles2)
## setDefaultNodeColor(cw, 'transparent')
## setNodeOpacityDirect(cw, names(pieData), 0)
## redraw(cw)


## ----desaturation--------------------------------------------------------
library(colorspace)
desaturatePercentage <- function(col, percentage=0.5){
  # this code is taken from the desaturate function in colorspace
  if (is.character(col) && (all(substr(col, 1L, 1L) == "#") & 
    all(nchar(col) %in% c(7L, 9L)))) {
    alpha <- substr(col, 8L, 9L)
    col <- substr(col, 1L, 7L)
    col <- hex2RGB(col)
    }
  else {
    col <- col2rgb(col, alpha = TRUE)
    alpha <- format(as.hexmode(col[4L, ]), width = 2L, upper.case = TRUE)
    alpha[alpha == "FF"] <- ""
    col <- RGB(t(col[1L:3L, ])/255)
  }
  col <- as(col, "polarLUV")
  col@coords[, 2L] <- col@coords[, 2L] * percentage
  col <- hex(col)
  return(col)
}

twoColors <- c("red", desaturatePercentage("red", 0.25))
pie(c(0.5, 0.5), col=twoColors)


## ----genMultiplePies-----------------------------------------------------
library(categoryComparePaperRev)

ccOpts <- new("ccOptions", listNames=c("T10_UP", "T10_DN", "T48_UP", "T48_DN", "T1_X", "T1_Y"), colorType="pie")

nTerm <- 1000

pieMatrix <- do.call(cbind, lapply(listNames(ccOpts), function(inName){
    baseCol <- rep(ccOpts@unsaturatedColor[inName], nTerm)
    baseCol[sample(nTerm, 500)] <- compareColors(ccOpts)[inName]
    baseCol
  }))

pieArea <- rep(1 / length(listNames(ccOpts)), length(listNames(ccOpts)))
names(pieArea) <- ""



