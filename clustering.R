## This script will do a Principal Component Analysis over a data set. After that kmeans algorithm will be implemented over the result.

library(raster)
library(rasterVis)
library(ncdf4)
library(maps)
library(maptools)

## Read the .nc data with raster

my_nc <- stack("my_nc.nc")


## 1. Firts function: do the PCA for the data we provide.

pca2kmeans <- function(x){
    data <- as.data.frame(x)
    datapca <-prcomp(data) ## PCA descomopsition
    ## I SELECT THE NUMBER OF MODES THAT I CONSIDER. 
    cumvar <- cumsum(datapca$sdev^2 / sum(datapca$sdev^2))
    b <- which(cumvar < 0.96 )
    c <- length(b)
    datapca2 <- data.frame(datapca$x[,1:c])
    return(datapca2)
}

## in case your data has NA -> data[is.na(data)] <- -999

datakm <- pca2kmeans(my_nc) ## The data frame resulting from the function is the one we are going to use for the kmeans.

## 2. 2nd function: do the kmeans algorithm n times. Every time will do it for 1:k clusters.

kmeansexp <- function(x, n, k){
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

resultado <- kmeansexp(datakm,n,k) ## you can especify n=number of running times and k=number of clusters. In my case I did n=500 and k=70

## resultado has a list with n elements. Yo need to select which one is the one you would like to chose. We use the CH index.

## 3. 3rd function to analyse the correct number of cluters.

CH <- function(x, k){
    CH_index <- c()
    for (i in 2:k) CH_index[i] <- ((nrow(x)-i)/(i-1))*(x[[i]]$tot.withinss)/(x[[i]]$betweenss)
return(CH_index)
                               }

Indice <- lapply(resultado, FUN='CH')
Indice <- do.call(cbind, Indice) 


## Once we have applied the index, we take the optimal partition and visualize it.

## How can I find the median vector on the distribution?

## Calc of the rmse of the lines like in the paper of Zagouras. (CALCULO DE LOS n RMSEc. Para ello, necesito calcular para cada exp n rmse y escoger el menor.)

rmse <- list()

ajuste <- function(x){
    
            for(j in 1:67){
        r1L <- lm(x[1:j+1]~c(1:j+1))
        r1R <- lm(x[j+2:70]~c(j+2:70))
    rmse[[j]] <- c(sqrt(sum((r1L$residuals)^2)), sqrt(sum((r1R$residuals)^2)))}
    return(rmse)
}

## I apply this 'ajuste' function to every element of 'Indice' (they are lists of kmeans results)

rmse_Exp <- apply(Indice, 2, FUN='ajuste')

## Ponderate the rmse from lm results.
ajustePond <- function(x){
    rmseT <- c()
        for(i in 1:67)
    rmseT[i] <-((i+1)-1)/((70)-1)*(x[[i]][1])+(70-(i+1))/(70-1)*(x[[i]][2])
    return(rmseT)
}


rmse_ExpP <- lapply(rmse_Exp, FUN='ajustePond')

## list of 'n' elements with smaller RMSE:

minRmse <- lapply(rmse_ExpP, FUN=function(x) unlist(x)) ## list of n elements with smaller RMSE

minRmse_w<- lapply(minRmse, FUN=function(x) min(x))

## Which of the n is on the middle

middle <- function(x) {
vectorMIN <- c(unlist(x))
mediano <- sort(vectorMIN)[50]

posicion <-which(vectorMIN == mediano) 
return(posicion)
}

middle(minRmse_w) ## Te da la posición de la partición mediana.


## Me voy a la posición 'middle(minRmse_w)' de la lista que contiene los RMSE para todos los ajustes. Tengo que ver qué porsición corresponde al valor mínimo:

kopt <- min(rmse_ExpP[['middle(minRmse_w)']]) ## en los corchetes tengo que poner la posición que he obtenido con middle(minRmse_w)
kpos <- which(rmse_ExpP[['middle(minRmse_w)']] == kopt) ## Lo mismo
kpos

## Una ve que tengo la partición óptima, tengo que representarla.

## El numero de clusters siempre es kpos+1

## partición ótima:

optima <- resultado[['middle']][['kpos']]

clusters <- raster(my_nc)
clusters <- setValues(clusters, optima$cluster)

levelplot(clusters)
