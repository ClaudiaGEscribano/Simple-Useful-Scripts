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
    datapca <- prcomp(data, center=TRUE, scale=TRUE) ## PCA

    ## I select the number of PC that I consider

    cumvar <- cumsum(datapca$sdev^2)/sum(datapca$sdev^2)
    b <- which(cumvar < 0.95)
    c <- length(b)
    datapca2 <- data.frame(datapca$x[,1:c])
    return(datapca2)
}

datahc <- pca2cluster(SISd1)

#############################################################
## 3. Hierarquical clustering
#############################################################

## I applied a hierarquical criterion in order to inizilitates the afterwards kmeans clustering.

dataMatrix <- as.matrix(datahc)
dataDist <- dist(dataMatrix) ## computes the distance matrix in order to use it in the hc algorithm

datahclus <- hclust(dataDist)

## I get the cut the clustering tree for k=2:70 clusters. I will decide the optimun number in next steps.

kut <- cutree(datahclus, k=2:70)

## I have to select the centroids of the clusters in order to use them in the kmeans clustering.
 
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

cluster.centroid <- function(datahc, data, n){ ## n es el numero de clusters max. menos uno.
     lapply(seq(from=1, to=n),
            FUN=function(x) lapply(seq(from=1, to=x+1),
                FUN=function(i){
                    ind  <- which(datahc[,x]==i) ## cambio aqui  i por x
                    r <- data[ind,]
                    mr <- rowMeans(r)
                    mr2 <- which(mr == sort(mr)[length(mr)/2])
                    where <- ind[mr2]
                    return(where)
                })
            )
 }

centroids <- cluster.centroid(kut, dataMatrix, 69)

a <- which(kut[,2] == 3)
b <- dataMatrix[a,]
c <- rowMeans(b)
d <- which(c == sort(c)[length(c)/2])
e <- a[d]
