## ----knitr, echo=FALSE, results="hide"-----------------------------------
library("knitr")
opts_chunk$set(tidy=FALSE,dev="png",fig.show="hide",
               fig.width=10,fig.height=8,
               message=FALSE)

## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"------------------
BiocStyle::latex()

## ----loadpackage, echo=FALSE---------------------------------------------
# in order to print version number below
library("flowMap")

## ----options, results='hide', echo=FALSE---------------------------------
options(digits=3, width=80, prompt=" ", continue=" ")

## ----dataPrepExample, echo=TRUE------------------------------------------
sam1 <- read.table(system.file("extdata/sample.txt",package="flowMap"),
                   header=T)
str(sam1)
table(sam1$id)

## ----displayData---------------------------------------------------------
sam1 <- read.table(system.file("extdata/sample.txt",package="flowMap"),
                   header=T)
sam2 <- read.table(system.file("extdata/sample.txt",package="flowMap"),
                   header=T)
table(sam1$id)
table(sam2$id)

## ----makeplotSelfruns1,dev='pdf',echo=TRUE,warning=FALSE-----------------
mat1 = sam1[sam1$id==1,]
mat2 = sam2[sam2$id==3,]

# combine events from the two cell populations to 
# make pooled data
mat = rbind(mat1,mat2)

# sample 100 events from the pooled data
sampleSize = 100

# among the 100 events, sample events from the two cell populations
# such that the ratio of the cell population membership is the same
# as that in the pooled data
nn1 = round(sampleSize*table(mat$id)[1]/nrow(mat))
nn2 = round(sampleSize*table(mat$id)[2]/nrow(mat))
submat = rbind(mat1[sample(nrow(mat1),nn1),],mat2[sample(nrow(mat2),nn2),])
colnames(submat)[5] = "sam"

# plot MST of the 100 events
g1 = makeFRMST(submat)
par(mar=c(0,0,0,0))
plot(g1$g,vertex.label.cex=0.01,
     layout=layout.fruchterman.reingold(g1$g))

## ----makeplotSelfruns2,dev='pdf',echo=TRUE,warning=FALSE-----------------
mat1 = sam1[sam1$id==4,]; mat2 = sam2[sam2$id==5,]
mat = rbind(mat1,mat2)
sampleSize = 100
nn1 = round(sampleSize*table(mat$id)[1]/nrow(mat))
nn2 = round(sampleSize*table(mat$id)[2]/nrow(mat))
submat = rbind(mat1[sample(nrow(mat1),nn1),],mat2[sample(nrow(mat2),nn2),])
colnames(submat)[5] = "sam"
g1 = makeFRMST(submat)
par(mar=c(0,0,0,0))
plot(g1$g,vertex.label.cex=0.01,layout=layout.fruchterman.reingold(g1$g))

## ----makeplotSelfruns3,dev='pdf',echo=TRUE,warning=FALSE-----------------
mat1 = sam1[sam1$id==6,]
mat2 = sam2[sam2$id==6,]
mat1$id=1; mat2$id=2
mat = rbind(mat1,mat2)
sampleSize = 100
nn1 = round(sampleSize*table(mat$id)[1]/nrow(mat))
nn2 = round(sampleSize*table(mat$id)[2]/nrow(mat))
submat = rbind(mat1[sample(nrow(mat1),nn1),],mat2[sample(nrow(mat2),nn2),])
colnames(submat)[5] = "sam"
g1 = makeFRMST(submat)
par(mar=c(0,0,0,0))
plot(g1$g,vertex.label.cex=0.01,layout=layout.fruchterman.reingold(g1$g))

## ----displayData2--------------------------------------------------------
sam1 <- read.table(system.file("extdata/sample.txt"
            ,package="flowMap"),header=T)
sam2 <- read.table(system.file("extdata/sample.txt"
            ,package="flowMap"),header=T)
table(sam1$id)
table(sam2$id)

## ----compareSampleSelf,dev='pdf',warning=FALSE---------------------------
res1 = getFRest(sam1,sam2,sampleMethod="proportional",sampleSize=100,
                ndraws=100,estStat="median",ncores=NULL)
res1@ww
library(gplots)
par(mar=c(0,0,0,0))
heatmapCols <- colorRampPalette(c("red","yellow","white","blue"))(50)
heatmap.2(res1@ww,trace="none",col=heatmapCols,symm=FALSE,dendrogram="none",
          Rowv=FALSE,Colv=FALSE,xlab="Sample 2",ylab="Sample 1")

## ----plotSelfpval,dev="pdf",echo=TRUE------------------------------------
library(gplots)
par(mar=c(0,0,0,0))
heatmapCols <- colorRampPalette(c("red","yellow","white","blue"))(50)
heatmap.2(res1@pNorm,trace="none",col=heatmapCols,symm=FALSE,dendrogram="none",
          Rowv=FALSE,Colv=FALSE,xlab="Sample 2",ylab="Sample 1")

## ----plotSelfpvalhist,dev="pdf",echo=TRUE--------------------------------
hist(res1@pNorm,xlab="log10 p-value histogram",main="")

## ----plotMultipval,dev="pdf",echo=TRUE-----------------------------------
resMulti = makeDistmat(samples=list(sam1,sam2),sampleSize=100,ndraws=100)
require(gplots)
par(mar=c(0,0,0,0))
heatmapCols <- colorRampPalette(c("red","yellow","white","blue"))(50)
heatmap.2(resMulti$distmat,trace="none",col=heatmapCols,symm=FALSE,dendrogram="none",
          Rowv=FALSE,Colv=FALSE)

## ----sessionInfo,results="asis",echo=FALSE, eval=TRUE--------------------
toLatex(sessionInfo())

## ----resetOptions, results='hide', echo=FALSE----------------------------
options(prompt="> ", continue="+ ")

## ----closeSockets,echo=FALSE---------------------------------------------
closeSockets <- function() {
    allCon <- showConnections()
    socketCon <- as.integer(rownames(allCon)[allCon[, "class"] == "sockconn"])
    sapply(socketCon, function(ii) close.connection(getConnection(ii)) )
}
closeSockets()

