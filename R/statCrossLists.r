#' summary statistic across elements of lists
#' 
#' Compute an estimate of F-R statistics across random draws of cell population comparison. 
# 'Central tendency measures such as median or mean are recommended.
#'  
#' @param obj an object of one or more lists. List elements must be matrices of sample dimensions
#' @param STAT the statistic to be computed for each entry in the matrix across lists. Default value: median.
#'     
#' @return est a matrix of the estimated statistic across random draws.
#'
#' @name statCrossLists
#' 
#' @export
statCrossLists <- function(obj,STAT) {
	# require(abind)
	combinedLists <- abind(obj,along=3)
	est <- apply(combinedLists,c(1,2),STAT)
	return(est)
}
