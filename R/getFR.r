#' FR test for a single cell population pair comparison
#'
#' Compute FR statistic for any two cell populations in flow cytometry data. Runtime of the FR test is 
#' quadratic in the number of events (nodes) totaling a single population pair comprison. 
#' 
#' @param xx1 events of a single cell population, organized in a matrix or a data.frame of events (rows) by 
#'          features (columns).
#' @param xx2 events of a single cell population, organized in a matrix or a data.frame of events (rows) by 
#'          features (columns).
#' 
#' @return ww FR statistic.
#' @return runs number of within-group subtrees (large number of runs indicate high degree of dissimilarity
#'              between the two cell populations being compared.
#' @return mu expected number of runs (when the two cell populations are similarly distributed).
#' @return sigma2 variance of runs.
#' @return pNorm p-values of the FR statitic assuming large sample asymptotic normal assumption.
#'
#' @examples
#' ## see vignettes
#'
#' @name getFR
#' 
#' @author Chiaowen Joyce Hsiao \email{joyce.hsiao1@@gmail.com}
#' 
#' @export
getFR <- function(xx1,xx2)
{
  xx <- rbind(xx1,xx2)
  distmat <- dist(xx, method = "euclidean",diag = T, upper=T)

  mstree <- ade4::neig2mat(ade4::mstree(distmat))
  mstree <- as.matrix(mstree)

  n1 <- dim(xx1)[1]
  n2 <- dim(xx2)[1]

  leftbottom <- mstree[(n1+1):(n1+n2),1:n1]
  rightup <- mstree[1:n1,(n1+1):(n1+n2)]
  runs <- sum(rightup)+1

  xsum <- colSums(mstree)

  C <- sum(xsum*(xsum-1))/2
  m <- n1
  n <- n2
  N <- n2+n1

  mu <- 2*m*n/N + 1
  sigma2 <- (2*m*n/(N*(N-1)))*((2*m*n-N)/N+(C-N+2)*(N*(N-1)-4*m*n+2)/((N-2)*(N-3)))
  ww <- (runs-mu)/sqrt(sigma2)
  pNorm <- pnorm(ww)

  list(ww=ww,runs=runs,mu=mu,sigma2=sigma2,pNorm=pNorm)
}

