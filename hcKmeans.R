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
 

cluster.centroid <- function(datahc, data, n){ ## n es el numero de clusters max. menos uno.
     lapply(seq(from=1, to=n),
            FUN=function(x) lapply(seq(from=1, to=x+1),
                FUN=function(i){
                    ind  <- which(datahc[,x]==i) ## cambio aqui  i por x
                    r <- data[ind,]
                    c <- colMeans(r)
                    return(c)
                })
            )
 }



centroids <- cluster.centroid(kut, dataMatrix, 69) ## contiene una lista de 70 elementos. Cada uno de los elementos es otra lista con los valores de de los centroides, que son vectores de dimensión d igual a las columnas de data.

a <- which(kut[,2] == 3)
b <- dataMatrix[a,]
c <- colMeans(b)
d <- which(c == sort(c)[length(c)/2])
e <- a[d]


############################################################################################
## 3. CENTERS MATRIX TO KMEANS
##########################################################################################


centros <- lapply(seq(from=1, to=69),
                  FUN=function(i) do.call(rbind, centroids[[i]]))


kmeansexp <- lapply(centros,
                    FUN=function(i) kmeans(dataMatrix, centers=i, iter.max=1000)) ## Me da un  warning.
    
## Puedo representar cualquiera de los experimentos de kmeanexp para explorar como son las particiones.

P <- raster(SISd)
P <- setValues(P, kmeansexp[[17]]$cluster)

levelplot(P)

## Lo que parece es que los centros que se dan ahora están más cerca de los buenos, por lo que el kmeans es mucho más rápido. Aún así, puede aparecer una no convergencia. LO que sería buneo es darle, después de estos centros iniciales, la opción de que pudiera inicializarse de nuevo.

###################################################################################################
## 4. OPTIMUN K
#################################################################################################

## In order to determine the optimun numnber of clusters, I can use several indexes.

library(clusterCrit)

load('data_from_PCA_30years.Rdata')
datakm_matrix <- as.matrix(datahc)

## Davies-Bouldin index

criterioDB <- lapply(seq(from=2, to=70),
                     FUN=function(x) intCriteria(datakm_matrix, kmeansexp[[x]]$cluster, 'Davies_Bouldin')
                     )

## Draw the index vs k

critDB_v <- as.vector(criterioDB)

k <- c(2:70)
matplot(k, critDB_v, type="l", main='DB index', xlab='k', ylab='Davies_Boudin')

##############################################################################################
## 5. REASIGNAR CELDAS
########################################################################################

## Tomo para hacer como ejemplo kmeansexp[[20]]$clusters

klus <- P
 
## utilizo el paquete raster con la funcion adjacent para calcular el valor de los vecinos.

neig <- matrix(c(1,1,1,
                 1,0,1,
                 1,1,1), ncol=3, byrow=TRUE)

ad <- adjacent(P, 1:ncell(P), directions=8, pairs=TRUE, id=TRUE, sorted=TRUE) ## directions 8 y neig es lo mismo

tb <- table(P[ad[,2]], P[ad[,3]])

## Una vez que obtengo la matriz ad tengo que comparar los valores de mi raster y de los vecinos.

## celdas que corresponden a los vecinos de la celda id:

## ad[ad[,1] == id, 3]

## Me da varios valores, todos los vecinos de la celda id. Salvo para los márgenes, debe haber 8 valores. ¿A qué cluster pertenecen estos vecinos?

## P[ad[ad[,1] == id, 3]]

## Para comparar con el valor del cluster de la celda que estoy analizando:

P[ad[ad[,1] == 1, 3]] != P[1] ## me dice las que son distintas

## me da valores TRUE o FALSE

## Lo que quiero saber es el valor de la celda id que tiene vecinos con valor distinto a ella

distintos <- which(P[ad[ad[,1] == 1, 3]] != P[1]) ## me dice cuantos de los vecinos son distintos. En este caso ninguno

if (distintos > length(ad[,1]==id)/2) {
    P[id] <- ## aquí tengo que asignar el valor del cluster que más se repita en los vecinos 
