library(raster)
library(rasterVis)

## Cargo los datos de mi nc. Contienen un año de datos de irradiancia diaria (W/m2).

SIS <- stack("SIS_complete_83_13.nc") ## objeto que contiene 11323 capas, una por día.

## Tomo algunas capas de este rasterStack para seleccionar solo un año:

SIS <- subset(SIS, 1:365) ## me quedo con un año de datos

## El paquete raster tienen un algebra sencilla:

## Multiplico por 24 para obtener la radiación (Wh/m2)

SIS <- SIS*24

## Para encontrar la media de la serie temporal en cada celda:

mediaPeriodo <- mean(SIS) ## contiene un raster con una única capa con la media de los 365 días

## Para encontrar la media de cada capa:

mediaEspacial <- cellStats(SIS, stat='mean') ## contiene 365 valores con la media de cada capa.

## Para aplicar una función "customizada" por celda se utiliza la función 'calc'

prueba <- calc(SIS, fun=function(x) {x*2+3}) ## A cada celda (y en cada paso de tiempo multiplica por 2 el valor y le suma 3. Esta función no tiene ningún significado pero es como prueba.
