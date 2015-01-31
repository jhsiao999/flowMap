#' Generate a similarity matrix quantifying similarity between cell populations across flow cytmetry samples
#'
#' Perform \code{\link{getFRest}} for multiple sample comparisons. For a two-sample comparison of n1 and n2 
#' cell populations, we generate a similarity matrix of dimension (n1+n2)-by-(n1+n2).
#'  
#' @param samples a list of data.frames or matrices of flow cytometry samples. Each data.frame or matrix
#'            consists of columns of features (flow cytometry channels), followed by a column of population
#'			 	membership (\emph{id}).
#' @param sampleMethod downsampling method, options include \emph{equalSize} or \emph{proportional}. Both methods 
#'                    sample events without replacement from the combined events in a single cell population pair 
#'                    comparison. Using \emph{equalSize}, each sample includes an equal number of events from the
#'                    two cell populations being compared. Using \emph{proportional}, the ratio of the event 
#'                    membership is same as the ratio of event membershiop prior to sampling. 
#'                    (default: \emph{proportional})
#' @param sampleSize specifies \emph{S},the number of events to be included in each sample. For \emph{equalSize} sampling 
#'                  , S/2 is sampled from each population. For \emph{proportional} sampling, the ratio of event membership
#'                  is the same as the ratio of event membership prior to sampling. (default: 200)
#' @param ndraws number of random samples. Runtime is linear in the number of random samples. (default: 200) 
#' @return distmat a matrix of estimated FR statistics across multiple samples. Rows and columns are indexed in the order
#'					in which each sample appears in the input list. For example, a two sample comparison entails rows and
#'					columns of \emph{distmat} indexed as 1.x follwed by 2.y where x is the Sample 1 population IDs, and
#'					2.y is the Sample 2 population IDs.
#'
#' @name makeDistmat
#' 
#' @author Chiaowen Joyce Hsiao \email{joyce.hsiao1@@gmail.com}
#' 
#' @export
makeDistmat <- function(samples,sampleMethod="proportional",sampleSize=200,ndraws=200) {
	
	require(Matrix)
	samples_norm <- samples

	nn <- length(samples_norm)

	distmatList <- lapply(1:nn,function(j) {

	  obj1 <- samples_norm[[j]]
	  npops1 <- unique(obj1$id)
	  nnrows <- length(npops1)

	  list0 <- lapply(1:nn,function(i) 
		{
		obj2 <- samples_norm[[i]]
		npops2 <- unique(obj2$id)
		nncolumns <- length(npops2)

		if (i < j) {
			res <- matrix(0,nrow=nncolumns,ncol=nnrows)
		} else if (i >= j) {
			res <- getFRest(XX1=obj2,
					XX2=obj1,sampleMethod=sampleMethod,
							sampleSize=200,estStat="median",ndraws)
#		    res <- -log10(res@pNorm)
		    res <- res@ww
		}
		rownames(res) <- paste(i,sort(npops2),sep=".")
		colnames(res) <- paste(j,sort(npops1),sep=".")
		return(res)
		}) 
	  list0 <- do.call(rbind,list0)
	  return(list0)
	})
	distmat <- do.call(cbind,distmatList)
	distmat <- forceSymmetric(distmat,uplo="L")

	list(distmat=as.matrix(distmat))
}

