## This script apply PCA and two clustering methods: first hierarquical and secondly kmeans.

library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(rgdal)
library(mapdata)
library(zoo)
library(ncdf4)

#################################################################
## 1. LOAD THE DATA
###############################################################

SISd <- stack("/home/claudia/Documentos/satelite/SIS_complete_83_13.nc")

## Tomo un año para probar

SISd1 <- subset(SISd, 1:365)

################################################################
## 2. PCA
###############################################################

pca2cluster <- function(x) {
    data <- as.data.frame(x)
    datapca <- prcomp(data) ## PCA

    ## I select the number of PC that I consider

    cumvar <- cumsum(datapca$sdev^2)/sum(datapca$sdev^2)
    b <- which(cumvar < 0.95)
    c <- length(b)
    datapca2 <- data.frame(datapca$x[,1:c])
    return(datapca)
}

datahc <- pca2cluster(SISd1)

#############################################################
## 3. Hierarquical clustering
#############################################################

## I applied a hierarquical criterion in order to inizilitates the afterwards kmeans clustering.

dataMatrix <- as.matrix(SISd1)
dataDist <- dist(dataMatrix) ## computes the distance matrix in order to use it in the hc algorithm

datahclus <- hclust(dataDist)

## I get the cut the clustering tree for k=2:70 clusters. I will decide the optimun number in next steps.

kut <- cutree(datahclus, k=2:70)

## I have to select the centroids of the clusters in order to use them in the kmeans clustering.
 
clust.centroid <- function(i, dataClus, data){
    apply(data, 2, FUN=function(x){
              ind <- (dataClus==i)
              centroids <- which(mean(x[ind]) == sort(mean(x[ind]))[nrow(x[ind])/2])
              return(centroids)
          }
          )
}

foo <- function(clusterdata, data, x){
    ind <- (which(clusterdata==x))
    centroids <- which(rowMeans(data[ind,])== sort(rowMeans(data[ind,]))[length(data[ind,]/2)])
    idx <- ind[centroids]
    return(idx)
}

z <- lapply(seq(from=1, to=6), FUN=function(x) foo(prueba,x))

a <- foo(prueba, dataMatrix,3)

prueba <- kut[,6]

C <- lapply(seq(1:6), FUN=function(x){
                ind <- which(prueba == x)
                r <- dataMatrix[ind,]
                mr <- rowMeans(r)
                mr2 <- which(mr == sort(mr)[length(mr)/2])
                where <- ind[mr2]}
            )

foo <- function(x){
    ind <- which(prueba==x)
    r <- dataMatrix[ind,]
    mr <- rowMeans(r)
    mr2 <- which(mr == sort(mr)[length(mr)/2])
    where <- ind[mr2]
    return(where)
}

C <- lapply(seq(1:6), FUN='foo')

D <- lapply(seq(from=1, to=69),
            FUN=function(i) lapply(seq(from=1, to=i+1),
                FUN=function(x){
                    ind  <- which(kut[,i]==x)
                    r <- dataMatrix[ind,]
                    mr <- rowMeans(r)
                    mr2 <- which(mr == sort(mr)[length(mr)/2])
                    where <- ind[mr2]
                    return(where)
                })
            )
