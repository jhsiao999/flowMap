#' FR tests to compare two flow cytometry samples 
#'
#' Perform FR tests to compare cell populations across two flow cytometry samples.
#' The FR statistics are estimated based on a downsampling scheme designed to optimize runtime as well as
#' precision and accuracy of the statistics. The downsampling scheme samples from the pooled data of 
#' events across a single cell population comparison. The events in the sample maintain the same cell 
#' population membership ratio as the events in the pooled data (\emph{proportional}) or include an
#' equal number of events from the two cell populations in the comparison.
#'  
#' @param XX1 a flow cytometry sample of cell populations, organized in a matrix or a data.frame of events (rows)
#'            by features (columns) where cell populaiton memberships are indexed in the last column by 
#'            a variable named id.
#' @param XX2 a flow cytometry sample of cell populations, organized in a matrix or a data.frame of events (rows)
#'            by features (columns) where cell populaiton memberships are indexed in the last column by 
#'            a variable named id.
#' @param sampleMethod downsampling method, options include \emph{equalSize} or \emph{proportional}. Both methods 
#'                    sample events without replacement from the combined events in a single cell population pair 
#'                    comparison. Using \emph{equalSize}, each sample includes an equal number of events from the
#'                    two cell populations being compared. Using \emph{proportional}, the ratio of the event 
#'                    membership is same as the ratio of event membershiop prior to sampling. 
#'                    (default: \emph{proportional})
#' @param sampleSize specifies \emph{S},the number of events to be included in each sample. For \emph{equalSize} sampling 
#'					, S/2 is sampled from each population. For \emph{proportional} sampling, the ratio of event membership
#'					is the same as the ratio of event membership prior to sampling.
#' @param i dummy variable to initialize parallel computing (\pkg{doParallel}).
#' 
#' @return wmat a list of matrices containing sample FR statistics for each XX1 by XX2 population comparisons.
#' @return runsmat a list of matrices containing sample runs for each XX1 by XX2 population comparisons.
#' @return mumat a list of matrices containing sample expected number of runs for each XX1 by XX2 population comparisons.
#' @return sigma2mat a list of matrices containing sample estimated variance of runs for each XX1 by XX2 population comparisons.
#' @return pNormat a list of matrices containing sample one-sided p-value associated with the FR statistics of 
#'                XX1 by XX2 population comparisons under the asymptotic normality assumption of the FR statistic.
#'
#' @name getFRmat
#' 
#' @author Chiaowen Joyce Hsiao \email{joyce.hsiao1@@gmail.com}
#' 
#' @export
getFRmat <- function(XX1,XX2,sampleMethod,sampleSize,i=NULL) { 

  	XX1 <- as.data.frame(XX1)
    XX2 <- as.data.frame(XX2)

    XX1.list <-split(XX1,f=XX1$id)
    XX2.list <-split(XX2,f=XX2$id)

    lvl1 <- length(XX1.list)
    lvl2 <- length(XX2.list)

    F1 <- dim(XX1.list[[1]])[2]-1
    F2 <- dim(XX2.list[[1]])[2]-1 

    XX1.list <- lapply(XX1.list,function(xx) {
      xx <- xx[,1:F1]
      return(xx)})

    XX2.list <- lapply(XX2.list,function(xx) {
      xx <- xx[,1:F2]
      return(xx)})


    wmat <- matrix(NA,lvl1,lvl2)
    colnames(wmat)=as.numeric(sort(unique(XX2$id)))
    rownames(wmat)=as.numeric(sort(unique(XX1$id))) 

    mumat <- wmat
    sigma2mat <- wmat
    runsmat <- wmat
    pNormat <- wmat

    # set up index array
    sps <- length(unique(XX1$id))
    ref.groups <- sort(unique(XX2$id))
    temp <- matrix(ref.groups,ncol=length(ref.groups),nrow=sps,byrow=T)
    iiTest <- data.frame(sp=sort(unique(XX1$id)),temp)

    testPops <- as.numeric(as.character(iiTest[,1]))

    for (jj in seq_along(testPops)) {

	  refPops <- iiTest[jj, c(2:ncol(iiTest)) ]

      for (kk in seq_along(refPops)) {

		j <- testPops[jj]
        k <- as.numeric(as.character( refPops )[kk])

        iirow <- which(rownames(wmat)==j)
        iicol <- which(colnames(wmat)==k)

        xx1 <- XX1.list[[iirow]]
        xx2 <- XX2.list[[iicol]]

        if (sampleMethod=="equalSize") {

		   groupSize <- sampleSize/2

		   if (nrow(xx1) > groupSize) {
           iisam1 <- sample(nrow(xx1),groupSize)
          } else {
		   iisam1 <- sample(nrow(xx1),groupSize,replace=TRUE)
          }
          if (nrow(xx2) > groupSize) {
           iisam2 <- sample(nrow(xx2),groupSize)
          } else {
           # iisam2 <- c(1:nrow(xx2))
           iisam2 <- sample(nrow(xx2),groupSize,replace=TRUE)
          }          
         mat1 <- xx1[iisam1,]
         mat2 <- xx2[iisam2,]
        } else  if (sampleMethod=="proportional") {
		   prop_temp <- round(nrow(xx1)/(nrow(xx1)+nrow(xx2)),2)
		   if (prop_temp==1) {prop_take=0.99}
		   if (prop_temp<1) {prop_take=prop_temp}
		   n1 <- sampleSize*prop_take
           n2 <- sampleSize-n1
           iisam1 <- sample(nrow(xx1),n1,replace=TRUE)
           iisam2 <- sample(nrow(xx2),n2,replace=TRUE) 
           mat1 <- xx1[iisam1,]
           mat2 <- xx2[iisam2,]
        } else if (sampleMethod=="none") {
            mat1 <- xx1
            mat2 <- xx2
        } else {
          message("warning: sampling method unspecified")
        }

        temp <- getFR(mat1,mat2)

        wmat[iirow,iicol] <- temp$ww
        mumat[iirow,iicol] <- temp$mu
        sigma2mat[iirow,iicol] <- temp$sigma2
        runsmat[iirow,iicol] <- temp$runs
        pNormat[iirow,iicol] <- temp$pNorm

      }      
    }

    list(wmat=wmat,runsmat=runsmat,mumat=mumat,sigma2mat=sigma2mat,pNormat=pNormat)  
}


