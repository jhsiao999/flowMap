closeSockets <- function() {
    allCon <- showConnections()
    socketCon <- as.integer(rownames(allCon)[allCon[, "class"] == "sockconn"])
    sapply(socketCon, function(ii) close.connection(getConnection(ii)) )
}
