## Este script es para calcular las correlaciones en las series de radiación diaria entre los clusters de la peninsula.

library(raster)
library(rasterVis)
library(maps)
library(ncdf4)
library(maptools)
library(mapdata)
library(rgdal)

load('mascarasClusters.Rdata') ## contiene una lista con las mascaras de los clusters

## medias por cluster

serieSIScluster <- function(x, k){ ## es el el raster de la peninsula con los datos de radiacion y k es el cluster que quiero analizar
    data <- mask(x, mascarasClusters[[k]])
    ave <- cellStats(data, stat='mean')
    return(ave)
}

rsds <- stack("x")

series <- list(1:6)
for (i in 1:6) series[[i]] <- serieSIScluster(rsds, i)

correlaciones <- function(x, a, b, p){
    a <- series[[a]]
    b <- series[[b]]
    ab <- cor.test(a,b,p)
    return(ab)
}

c1c2 <- correlaciones(rsds, 1, 2, 0.95)

## Agrupo por meses:

## Enero

monthIndex <- function(k){            
c <- c(1:20)
j <- k
for (i in 1:20){
    c[i] <- j
    j <- j+12}                         
return(c)
}
            
lista <- list()
for (i in 1:12) lista[[i]] <- monthIndex(i)

## tengo que hacer la media por meses

library(zoo)

idx <- seq(as.Date("1989-01-01"), as.Date("2008-12-31"), 'day')
x <- setZ(x, idx)

SISmm240 <- zApply(x, by=as.yearmon, fun='mean')

Enero <- subset(SISmm240, lista[[1]])

series <- list(1:6)
for (i in 1:6) series[[i]] <- serieSIScluster(Enero, i)

c1c2 <- correlaciones(Enero, 1, 2, 0.95)

## Se hace la misma operación para las correlaciones entre todos los clusters y los meses que quiera.
