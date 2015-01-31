#' Plot the mimumum spanning tree of a cell population comparison
#'
#' Visualize minimum spanning tree (MST)  of the pooled data combined from the two cell populations. MST is the basis of 
#' the FR statistics. Two populations are similar if their respective events congregate with events of the same population
#' membership. Runs is calculated as the number of edges connecting nodes of different cell population membership plus 1 (
#' (or equivalently, the number of subtrees of homogeneous cell population membership).
#' 
#' @param mat a matrix or data.frame of the pooled data of a cell population pair. The rows contain events in the pooled data.
#'            The columns are the expression markers. The population membership is indexed in the last column of the matrix or 
#'			  data.frame. The membership IDs are for plotting purposes, and hence need to uniquely identify the two cell
#'            populations in the pooled data.
#' @param node.colors colors to label cell population membership of the events. (default: c("blue","red"))
#' 
#' @return g an \emph{igraph} object that contains graph adjacency matrix of the minimum spanning tree calculated from
#'				euclidean distance between nodes in the pooled data combined from the two cell populations.
#' @return gall an \emph{igraph} object that contains graph adjacency matrix containing euclidean distance between
#'				nodes in the pooled data from the two cell population. Complete graph is assumed here.
#' @return ww FR statistics for the cell population comparison (formula: (runs-mu)/sqrt(sigma2))
#' @return runs observed number of runs for the cell population comparison
#' @return mu expected number of runs for the cell population comparison
#' @return sigma2 variance of runs for the cell population comparison
#' @return pNorm p-value of the FR statistic for the cell population comparison
#' @return distmat a matrix defining Euclidean distance between events (nodes) across the two cell populations 
#' @return mstree minimum spannin tree defined in an adjacency matrix. 1 indicates an edge between nodes and 0 indicates 
#'					non-existent edge.
#' @return C number of edge pairs (pairs of edges that share a common node) in the minumum spanning tree
#' @return m number of events in the first cell population
#' @return n number of events in the second cell population
#'
#' @examples
#' ## see vignettes
#'
#' @name getFRest
#' 
#' @author Chiaowen Joyce Hsiao \email{joyce.hsiao1@@gmail.com}
#' 
#' @export
makeFRMST <- function(mat,node.colors=c("blue","red")) {

	Nevents = table(mat$sam)

	distmat <- as.matrix(dist(mat[,-ncol(mat)]),upper=T,diag=T)

	mstree <- ade4::neig2mat(ade4::mstree(as.dist(distmat)))
	mstree <- as.matrix(mstree)

	library(igraph)

	sample.labels <- c(rep(1,Nevents[1]),rep(2,Nevents[2]))

	g <- graph.adjacency(distmat*mstree,mode="undirected",weighted=TRUE,diag=FALSE,
	add.colnames=TRUE,add.rownames=TRUE)
	set.edge.attribute(g,name="weight",value=E(g)$weight^15)
	E(g)$color <- "black"
	E(g)$width <- 3
	V(g)$color <- node.colors[as.numeric(sample.labels)]
	V(g)$size <- 5

	gall <- graph.adjacency(distmat,mode="undirected",weighted=TRUE,diag=FALSE,
	add.colnames=TRUE,add.rownames=TRUE)
	set.edge.attribute(gall,name="weight",value=E(gall)$weight^15)
	E(gall)$color <- "grey"
	V(gall)$size <- 4
	V(gall)$color <- node.colors[as.numeric(sample.labels)]

	# compute statistics
	xx1 <- mat[mat$sam==1,]; xx2 <- mat[mat$sam==2,]
	n1 <- nrow(xx1)
	n2 <- nrow(xx2)

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

	list(g=g,gall=gall,ww=ww,runs=runs,mu=mu,sigma2=sigma2,pNorm=pNorm,
	distmat=distmat,mstree=mstree,C=C,m=m,n=n)
}





