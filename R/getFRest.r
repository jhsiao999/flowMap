#' FR tests to compare two flow cytometry samples 
#'
#' Estimate FR statistics comparing cell populations across two flow cytometry samples. The estimates are 
#' generated from \code{\link{getFRmat}} and are dependent on the number of random draws as well as the size
#' of each random sample. For every cell population pair comparison, we take the median of the FR statistics
#' across all random samples as the estimated similarity. Parallel computing is used to minimize the runtime
#' of the algorithm (\pkg{doParallel}). Users can specify the number of cores to be used depending on the 
#' available computing power. Default uses all available processing cores in the system. 
#' 
#' @param XX1 a flow cytometry sample of cell populations, organized in a matrix or a data.frame of events (rows)
#'            by features (columns) where cell populaiton memberships are indexed in the last column by 
#'            a variable named id.
#' @param XX2 a flow cytometry sample of cell populations, organized in a matrix or a data.frame of events (rows)
#'            by features (columns) where cell populaiton memberships are indexed in the last column by 
#'            a variable named id.
#' @param ndraws number of random samples. Runtime is linear in the number of random samples. (default: 200) 
#' @param sampleMethod downsampling method, options include \emph{equalSize} or \emph{proportional}. Both methods 
#'                    sample events without replacement from the combined events in a single cell population pair 
#'                    comparison. Using \emph{equalSize}, each sample includes an equal number of events from the
#'                    two cell populations being compared. Using \emph{proportional}, the ratio of the event 
#'                    membership is same as the ratio of event membershiop prior to sampling. 
#'                    (default: \emph{proportional})
#' @param sampleSize specifies \emph{S},the number of events to be included in each sample. For \emph{equalSize} sampling 
#'                  , S/2 is sampled from each population. For \emph{proportional} sampling, the ratio of event membership
#'                  is the same as the ratio of event membership prior to sampling. (default: 200)
#' @param estStat statistic that used to estimate FR statistic of each population pair comparison across random samples.
#'                  (default: median)
#' 
#' @return wmat a matrix of estimated FR statistics for each XX1 by XX2 population comparisons.
#' @return runsmat a matrix of estimated runs for each XX1 by XX2 population comparisons.
#' @return mumat a matrix of estimated expected number of runs for each XX1 by XX2 population comparisons.
#' @return sigma2mat a matrix of estimated variance of runs for each XX1 by XX2 population comparisons.
#' @return pNormat a matrix of one-sided p-values of the estimated FR statistic for each XX1 by XX2 population
#'                comparisons under the asymptotic normality assumption of the FR statistic.
#'
#' @examples
#' ## see vignettes
#'
#' @name getFRest
#' 
#' @author Chiaowen Joyce Hsiao \email{joyce.hsiao1@@gmail.com}
#' 
#' @export
getFRest <- function(XX1,XX2,sampleMethod="proportional",
		sampleSize=200,estStat="median",ndraws=200,ncores=NULL) {

	if (is.null(ncores)) {
	  ncores <- detectCores()
	  registerDoParallel(cores=ncores)
	  message(paste("computing on",ncores,"cores"))
	} else {
	  registerDoParallel(cores=ncores)
	  message(paste("computing on",ncores,"cores"))
	}
	
    message("computing FR statistics between sample ... \n")
	mat <- foreach(i=1:ndraws) %dopar% getFRmat(XX1,XX2,sampleMethod=sampleMethod,
        sampleSize=sampleSize,i)
    closeSockets()

    npop1 <- length(unique(XX1$id))
    npop2 <- length(unique(XX2$id))


    wmat <- lapply(mat,"[[",1) 
    wmat <- statCrossLists(wmat,estStat)
    colnames(wmat)=as.numeric(as.character(sort(unique(XX2$id))))
    rownames(wmat)=as.numeric(as.character(sort(unique(XX1$id))))

    runsmat <- lapply(mat,"[[",2) 
    runsmat <- statCrossLists(runsmat,estStat)

    mumat <- lapply(mat,"[[",3) 
    mumat <- statCrossLists(mumat,estStat)

    sigma2mat <- lapply(mat,"[[",4) 
    sigma2mat <- statCrossLists(sigma2mat,estStat)

    pNormat <- lapply(mat,"[[",5) 
    pNormat <- statCrossLists(pNormat,estStat)

    new("FRstats",
    	XX1=XX1,
    	XX2=XX2,
		sampleMethod=sampleMethod,
		sampleSize=sampleSize,
    	ncores=ncores,
        ndraws=ndraws,
    	npop1=npop1,
    	npop2=npop2,
    	pop1Labels=sort(unique(XX1$id)),
    	pop2Labels=sort(unique(XX2$id)),
    	ww=wmat,
    	runs=runsmat,
    	mu=mumat,
    	sigma2=sigma2mat,
    	pNorm=pNormat)
}





