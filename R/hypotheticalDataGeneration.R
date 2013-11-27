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