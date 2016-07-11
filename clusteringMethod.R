## This script will calculate kmeans algorithm to do a clustering partition over a dataset.

library(raster)
library(rasterVis)
library(maps)
library(maptools)
library(mapdata)
library(rgdal)
library(zoo)
library(ncdf4)

########################################################################################
## 1. LOAD DE DATA
########################################################################################

## 30 years of daily irradiation data from CMSAF

SISd <- stack("/home/claudia/Documentos/satelite/SIS_complete_83_13.nc")

## creo el objeto de polígonos que representa la península:

ext <- as.vector(extent(SISd))
boundaries <- map('worldHires', region=c('Spain','Portugal'),fill=TRUE, xlim=ext[1:2], ylim= ext[3:4], plot=FALSE)
boundaries$names
IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
boundaries_sp<- map2SpatialPolygons(boundaries, IDs=IDs, proj4string=CRS(projection(SISd)))
 
linea <- as(boundaries_sp, "SpatialLines") ## just the coast line of the IP

#########################################################################################
## 2. KMEANS k=1:k y algoritmo n veces
#########################################################################################


## 1. First function: do the PCA for the data we provided.

pca2kmeans <- function(x){
    data <- as.data.frame(x)
    datapca <-prcomp(data) ## PCA descomopsition
    ## I SELECT THE NUMBER OF MODES THAT I CONSIDER. 
    cumvar <- cumsum(datapca$sdev^2 / sum(datapca$sdev^2))
    b <- which(cumvar < 0.95 )
    c <- length(b)
    datapca2 <- data.frame(datapca$x[,1:c])
    return(datapca2)
}

datakm <- pca2kmeans(SISd)

##  2nd function: do the kmeans algorithm n times. Every time will do it for 1:k clusters.

kmeansexp <- function(x, n, k){ #datakm, n número de veces que repito el algoritmo, k número de clusters.
    expkm <- list(1:n)
    
     foo <- function(x, k){
        km <- list(1:k)

        for (i in 2:k){
        km[[i]] <- kmeans(x, i, nstart=25, iter.max=2000)}
        return(km)
    }
     for (i in 1:n){
        expkm[[i]] <- foo(x, k)} 
    return(expkm)
}

resultado <- kmeansexp(datakm,500,70) ## you can especify n=number of running times and k=number of clusters. In my case I did n=500 and k=70


############################################################################################
## 3. CH validation index and SSE elbow criterion
###########################################################################################

## To validate the results and to find the k=optimun I compute a validation index -> ch

##índice CH

CH <- function(x){
    CH_index<- c()
    for (i in 2:70) CH_index[i] <- ((14352-i)/(i-1))*(x[[i]]$betweenss)/(x[[i]]$tot.withinss)
    return(CH_index)
}

## sse (Normally we can find an elbow in sse, but it is not easy with radiation data). To compare with CH index we compute also the elbow criterion.

SSE <- function(x){
twss <- c()
for (i in 2:70) twss[i] <- x[[i]]$tot.withinss
return(twss)
}


SSE2 <- function(x){
    twss <- c()
    for (i in 2:70) twss[i] <- x[[i]]$tot.withinss

    tss <- c()
    for(i in 2:70) tss[i] <- x[[i]]$totss
    
    vari <- 1-twss/tss
    return(vari)
}
    ## para saber el número de clusters mínimo para explicar un porcentaje de variacion añadimos esta parte de la función.
    
    k <- min(which(vari>0.7))
    return(k)
}


##########################################################################################

## Apply the selected index (CH or SSE/SSE2) to the results of kmeans 500 runs experiments.

criterio <- lapply(resultado, FUN='SSE') ## calculo el índice para cada experimento.
criterio <- do.call(cbind, criterio) ## put in columns the 'x' experiments. There are 70 rows with  tot.withinss ## matriz


## We can represent the variance in data explained depending on the number of cluster selected. Also it is good to show how the CH index varies with the number of clusters.

k <- c(1:70)
matplot(k, criterio, type='l', main='Variance explained depending oin the number of clusters')
abline(v=min(which(criterio>0.7))) ## this line show the number of clusters needed to explain a percentage of variance in data.

criterioCH <- lapply(resultado, FUN='CH')
criterioCH <- do.call(cbind, criterioCH)

matplot(k, criterioCH, type='l', main='CH index') 

## We make use of the clusterCrit package to compute more indexes. In addition, L-method is applied with 'ajuste' and 'ajusteExp' to decide the optimun number of clusters.

library(clusterCrit)

load('data_from_PCA_30years.Rdata')
datakm_matrix <- as.matrix(datakm)


criterioDB <- lapply(resultado, FUN=function(x)
    lapply(seq(from=2, to=70), FUN= function(i)
        intCriteria(datakm_matrix, x[[i]]$cluster, 'Davies_Bouldin')
           )
)

load('DB_index_example.Rdata')

k <- c(2:70)
matplot(k, crit_v, type="l", main='DB index', xlab='k', ylab='RMSEt')

ajuste <- ajuste(crit_v)
ajuste_t <- ajustePond(ajuste)

matplot(k, ajuste_t, type="l",main="DB index min for RMSEt", xlab='k', ylab='RMSEt') ##26
abline(v=which(ajuste_t == min(ajuste_t)), col=2) 


## ajuste para SSE

criterioExample <- criterio[,1]

ajusteSSE <- ajuste(criterioExample)
ajusteSSE_t <- ajustePond(ajusteSSE) ## 17 (al que habrá que sumar 2, 19)
 
k <- c(1:67)
matplot(k, ajusteSSE_t, type="l",main="SSE index min for RMSEt", xlab='k', ylab='RMSEt') ##26
abline(v=which(ajusteSSE_t == min(ajusteSSE_t)), col=2) 

## Probamos otro índice, el C por ejemplo

prueba <- resultado[[1]]
critC <- lapply(seq(from=2, to=70), FUN=function(i)
    intCriteria(datakm_matrix, prueba[[i]]$cluster, 'C_index')
)

criterioC <- lapply(resultado, FUN=function(x)
    lapply(seq(from=2, to=70), FUN= function(i)
        intCriteria(datakm_matrix, x[[i]]$cluster, 'C_index')
           )
)



critC_v <- as.vector(unlist(critC))
ajusteC <- ajuste(critC_v)
ajusteC_t <- ajustePond(ajusteC) ## 18 (al que habrá que sumar 2, 20)
 

k <- c(1:67)
matplot(k, ajusteC_t, type="l",main="C index min for RMSEt", xlab='k', ylab='RMSEt') 
abline(v=which(ajusteC_t == min(ajusteC_t)), col=2)

k <- c(2:70)
matplot(k, critC_v, type="l",main="C index vs k", xlab='k', ylab='RMSEt') 

## Probamos el índice Sihoulette

prueba <- resultado[[1]]
critS <- lapply(seq(from=2, to=70), FUN=function(i)
    intCriteria(datakm_matrix, prueba[[i]]$cluster, 'Silhouette')
)

critS_v <- as.vector(unlist(critS))
ajusteS <- ajuste(critS_v)
ajusteS_t <- ajustePond(ajusteS) 
 

k <- c(1:67)
matplot(k, ajusteS_t, type="l",main="Sihoulette index min for RMSEt", xlab='k', ylab='RMSEt') 
abline(v=which(ajusteS_t == min(ajusteS_t)), col=2) ## (Al que habrá que sumar 2, 16)

k <- c(2:70)
matplot(k, critS_v, type="l",main="Sihoulette index vs k", xlab='k', ylab='RMSEt')

## índice Dunn

prueba <- resultado[[1]]
critDu <- lapply(seq(from=2, to=70), FUN=function(i)
    intCriteria(datakm_matrix, prueba[[i]]$cluster, 'Dunn')
)

critDu_v <- as.vector(unlist(critDu))
ajusteDu <- ajuste(critDu_v)
ajusteDu_t <- ajustePond(ajusteDu) 
 

k <- c(1:67)
matplot(k, ajusteDu_t, type="l",main="Dunn index min for RMSEt", xlab='k', ylab='RMSEt') 
abline(v=which(ajusteDu_t == min(ajusteS_t)), col=2) ## (Al que habrá que sumar 2, 16)

k <- c(2:70)
matplot(k, critDu_v, type="l",main="Dunn index vs k", xlab='k', ylab='RMSEt')



##############################################################
## 4. OPTIMUN PARTITION
###############################################################

## We have to decide which is the optimum partition based on a criteria and making use of the validity indexes.

## Following the methodology of the Zagouras' paper, after having calculated the CH index, we will do a linear fitting that will help us to find the optimun number of clusters.

## In the graph CH vs k, we will fit the data to 2 lines: one for the left side and one for the right side. This two lines will have at least 2 points. The left side line will increase its number of point each time and the rigth side line will decrease it. We calculate the RMSE for the 2 lines at each fitting and after that we calculate a total error considering the two RMSE. The minimum of the RMSEtotal give us the optimun number of clusters.

## We do that for the 500 clustering experiments.

rmse <- list()

ajuste <- function(x){ ## ajusta a 2 rectas los puntos de la gráfica CH vs k.
            for(j in 1:67){
        r1L <- lm(x[1:j+1]~c(1:j+1))
        r1R <- lm(x[j+2:70]~c(j+2:70))
    rmse[[j]] <- c(sqrt(sum((r1L$residuals)^2)), sqrt(sum((r1R$residuals)^2)))}
    return(rmse)
}


## I apply this 'ajuste' function' to every element of 'elbow' (they are lists of kmeans results)

rmse_Exp <- apply(criterioCH, 2, FUN='ajuste')

## Ponderate the rmse from lm results.

ajustePond <- function(x){ ## Calculating the total error.
    rmseT <- c()
        for(i in 1:67)
    rmseT[i] <-((i+1)-1)/((70)-1)*(x[[i]][1])+(70-(i+1))/(70-1)*(x[[i]][2])
    return(rmseT)
}


rmse_ExpP <- lapply(rmse_Exp, FUN='ajustePond')

## In order to visualize how the RMSEtotal is evolving with the number of clusters, we can select one of the 500 experiments as an example and represent the graph: 

k <- c(3:69) ## First point in the graph corresponf to 3 clusters, so to the minimum we have to add 2 points. (Theres is no CH for 1 cluster and no RMSEtfor 2, it start in 3).

minimo <- min(rmse_ExpP[[24]])
kop <- which(rmse_ExpP[[24]]== minimo) 
plot(k, rmse_ExpP[[24]], type='p', xlab='k=number of clusters', ylab='RMSEt', main='Total RMSE vs number of clusters')
abline(v=19,col=2) 

crit <- lapply(prueba, FUN=function(x)
    intCriteria(data,x$cluster, 'Davies_Bouldin'))

crit <- c()
for (i in 2:70) crit[i] <- intCriteria(data, prueba[[i]]$cluster,'Calinski_Harabasz')
    
## Once we have calculated the total error, we can find which of the 500 kmeans experiments will be the selected. We have decided that the median partition from the 500 experiments will be selected.


## We can see how it is the rmse for one of the experiments, the 24.

## list of 500 elements with smaller RMSE:

minRmse <- lapply(rmse_ExpP, FUN=function(x) unlist(x)) ## list of 500 elements with smaller RMSE
minRmse_w<- lapply(minRmse, FUN=function(x) min(x))

## Which of the 500 is on the middle

middle <- function(x) {
vectorMIN <- c(unlist(x))
mediano <- sort(vectorMIN)[50]

posicion <-which(vectorMIN == mediano) 
return(posicion)
}

middle(minRmse_w) ## 212 es el experimento en la mediana

## Me voy a la posición 212  de la lista que contiene los RMSE para todos los ajustes. Tengo que ver que posición corresponde al valor minimo y esto me dará k para est aparticion mediana.

## NUMBER OF OPTIMAL CLUSTERS IN THE MEDIAN PARTITION
  
medianPtt <- rmse_ExpP[[86]]
kopt <- min(rmse_ExpP[[86]])
kpos <- which(rmse_ExpP[[86]] == kopt)
kpos ## 5

 
kopt <- min(rmse_ExpP[[327]])
kpos <- which(rmse_ExpP[[327]] == kopt)

       
######################################

## Según hemos visto en el índice CH, la partición mediana es la numero 212 y 5 (6) es el numero de clusters optimo.

## cuando hago CH multiplicado por un resultado, obtengo q el numero optimo de clusters es 6 

## Cuando utilizamos para la validacion SSE, obtenemos que la particion mediana es la 19 y que el número óptimo de clusters es 17.

## OJO! CUANDO UTILIZO BIEN EL ÍNDICE CH ME SALEN 5 CLUSTERS

#############################################################################################################################
## 4. TOMO LA PARTICIÓN ÓPTIMA
###########################################################################################################################

optimK <- resultado[[86]][[17]]

## TOMO UNA PARTICION Y LA CONVIERTO EN RASTER

cluster6sat <- raster(SISd)
cluster6sat <- setValues(cluster6sat, optimK$cluster)

levelplot(cluster6sat, margin=FALSE, main=list('CM-SAF cluster partition k=8', cex=2))+layer(sp.lines(linea))

dev.copy(png, file='optimalpartitionCMSAF_30_8c.png', width=800, height=800)
dev.off()

###################################################################################################
## 5. Converting every cluster in a SpatialPolygons dataframe
##################################################################################################

Polygons <- list()
for (i in 1:6) Polygons[[i]] <- rasterToPolygons(cluster6sat, fun=function(x) {x==i})
levelplot(mask(cluster6sat, Polygons[[1]]), margin=FALSE, main='1')+layer(sp.lines(linea))
levelplot(mask(cluster6, Polygons[[2]]), margin=FALSE, main='2')+layer(sp.lines(linea))
levelplot(mask(cluster6, Polygons[[3]]), margin=FALSE, main='3')+layer(sp.lines(linea))
levelplot(mask(cluster6, Polygons[[4]]), margin=FALSE, main='4')+layer(sp.lines(linea))
levelplot(mask(cluster6, Polygons[[5]]), margin=FALSE, main='5')+layer(sp.lines(linea))
levelplot(mask(cluster6, Polygons[[6]]), margin=FALSE, main='6')+layer(sp.lines(linea)) 


ks <- list(1:6)
for (i in 1:6) ks[[i]] <- mask(cluster6sat, Polygons[[i]]) 
ksB <- lapply(ks, FUN=function(x) mask(x, boundaries_sp)) ## lista que contiene las mascaras de los clusters dentro de la peninsula

## guardo la máscara de clusters del satélite

save(ksB, file='mascaraClustersSat.Rdata')



################### 13 abril  2016 ######################################################################



## Tomo los datos de radiacion para hacer los calculos que necesito. Para ello, tengo que transformar en .nc SISS que es el rasterStack corregido y así poder utilizar las cdo.

## Al rasterStack SISS es necesario darle nombres primero

idx <- seq(as.Date('1989-01-01'), as.Date('2008-12-30'), 'day')
names(SISS) <- idx
SISS <- setZ(SISS, idx)

years <- function(x) as.numeric(format(x,'%y'))
SISym <- zApply(SISS, by=years, fun=mean)
 

## INTERANUAL STATISTICAL ANALYSIS ##

SISym <- stack("rsdsym.nc", varname="rsds")
SISym <- SISym*8760/1000
idx <- seq(as.Date('1989-01-01'), as.Date('2008-12-31'), 'year')
SISym <- setZ(SISym, idx)
names(SISym) <- idx

SISymd <- stack("rsdsym.nc", varname="rsds")
idx <- seq(as.Date('1989-01-01'), as.Date('2008-12-30'), 'year')
SISymd <- setZ(SISymd, idx)
names(SISymd) <- idx


## máscara de clusters ksB sobre la radiacion media anual

SISym_clusters <- list(1:6)
for (i in 1:6) SISym_clusters[[i]] <- mask(SISym, ksB[[i]])

## CV cluster

years <- function(x) as.numeric(format(x,'%y'))

VariabilidadInteranual <- function(x) {
    sd <- calc(x, fun=function(x) sd(x))
    media <- mean(x, na.rm=TRUE)
    COV <- sd/media
    return(COV)
}

CV <- lapply(SISym_clusters, FUN='VariabilidadInteranual') 

## Hago algunas operaciones más para cada cluster:

UNO <- CV[[1]]

levelplot(mask(UNO, boundaries_sp))+layer(sp.lines(linea))
xyplot(UNO)

#######################################################

mediaCV <- lapply(CV, FUN=function(x) cellStats(x, stat='mean'))
 
CV_layers <- stack(CV[[1]], CV[[2]], CV[[3]], CV[[4]], CV[[5]], CV[[6]])
names(CV_layers) <- c("cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6")

boxplot(CV_layers, range=70, main='media CV interanual', col=2)
bwplot(CV_layers)
densityplot(CV_layers)

## Hay diferencias entre los clusters del modelo corregido y los del satelite.

## MONTHLY STATISTICAL ANALYSIS ##

library(zoo)

##240 meses
SISmm240 <- zApply(SISS, by=as.yearmon, fun='mean') ##raster with the 240 monthly daily global radiation mean.


##12 meses
month <- function(x) as.numeric(format(x, '%m'))

SISmm <- zApply(SISS, by=month, fun='mean')
names(SISmm) <- month.abb


idx2 <- seq(as.Date('1989-01-01'), as.Date('2008-12-31'), 'month')

## Pasamos la máscara de clusters.

SISmm240_clusters <- list(1:6)
for (i in 1:6) SISmm240_clusters[[i]] <- mask(SISmm240, ksB[[i]])
for (i in 1:6) SISmm240_clusters[[i]] <- setZ(SISmm240_clusters[[i]], idx2)


VariabilidadMensual <- function(x) {
    sd <- zApply(x, by=month, fun='sd')
    media <- zApply(x, by=month, fun='mean')
    COV <- sd/media
    return(COV)
}

CV_mon <- lapply(SISmm240_clusters, FUN='VariabilidadMensual')

mediaCV_mon<- lapply(CV_mon, FUN=function(x) cellStats(x, stat='mean'))

zz <- lapply(mediaCV_mon, FUN=function(x) zoo(x))
zzz <- do.call(cbind, zz)
colnames(zzz) <- c("1", "2", "3", "4", "5", "6")

xyplot(zzz, superpose=TRUE, auto.key=list(space='right'), main='Monthly CV 6 clusters')

CV_mon_stack <- stack(CV_mon[[1]],CV_mon[[2]],CV_mon[[3]],CV_mon[[4]],CV_mon[[5]],CV_mon[[6]])
levelplot(CV_mon_stack, layers=1:12)+layer(sp.lines(linea))
boxplot(CV_mon_stack[[24:36]])

## DAILY STATISTICAL ANALYSIS##

## Cálculo de la variabilidad diaria.

idx <- seq(as.Date('1989-01-01'), as.Date('2008-12-30'), 'day')
SISS <- setZ(SISS, idx)


## Esto no está bien, hay que corregirlo porque los días los coge del 1 al 31. No toma un numero distinto para cada dia del año. VOy a hacer este cálculo con las cdo.

day <- function(x) as.numeric(format(x, '%d'))


VariabilidadDiaria <- function(x) {
    sd <- zApply(x, by=day, fun='sd')
    media <- zApply(x, by=day, fun='mean')
    COV <- sd/media
    return(COV)
}

SISd_clusters <- list(1:6)
for (i in 1:6) SISd_clusters[[i]] <- mask(SISS, ksB[[i]])
for (i in 1:6) SISd_clusters[[i]] <- setZ(SISd_clusters[[i]], idx)

CV_day <- lapply(SISd_clusters, FUN='VariabilidadDiaria')

mediaCV_day<- lapply(CV_day, FUN=function(x) cellStats(x, stat='mean'))

### Haciendo el calculo con las cdo. 

std_daily <- stack("std_dailyCMSAF.nc", varname='SIS')
mean_daily <- stack("mean_dailyCMSAF.nc", varname='SIS')

CV_daily <- std_daily/mean_daily

SISd_clusters <- list(1:6)
for (i in 1:6) SISd_clusters[[i]] <- mask(CV_daily, ksB[[i]])
idx3 <- seq(as.Date('2008-01-01'), as.Date('2008-12-31'), 'day')

for (i in 1:6) SISd_clusters[[i]] <- setZ(SISd_clusters[[i]], idx3)

## media de cada cluster:

mediaCV_day<- lapply(SISd_clusters, FUN=function(x) cellStats(x, stat='mean')) ## obtengo 6 series temporales de 366 días


zz <- lapply(mediaCV_day, FUN=function(x) zoo(x))
zzz <- do.call(cbind, zz)
colnames(zzz) <- c("1", "2", "3", "4", "5", "6")

xyplot(zzz, superpose=TRUE, auto.key=list(space='right'), main='Daily CV 6 clusters')
dev.copy(bmp, file='dailyCV_6CLUSTERS_cmsaf.bmp', width=800, height=600)
dev.off()

xyplot(zzz, superpose=FALSE)
dev.copy(bmp, file='dailyCV_6CLUSTERS_500_2.bmp', width=800, height=600)
dev.off()

summary(zz[[1]])

xyplot(zz, superpose=TRUE)
lapply(mediaCV_day, FUN=function(x) histogram(x))

##########################################################################################

## REPRESENTACIONES DE CUADRANTES


