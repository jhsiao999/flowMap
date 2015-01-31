### R code from vignette source 'flowMap.Rnw'

###################################################
### code chunk number 1: options (eval = FALSE)
###################################################
## options(width=70)


###################################################
### code chunk number 2: closeSockets
###################################################
closeSockets <- function() {
    allCon <- showConnections()
    socketCon <- as.integer(rownames(allCon)[allCon[, "class"] == "sockconn"])
    sapply(socketCon, function(ii) close.connection(getConnection(ii)) )
}
closeSockets()


