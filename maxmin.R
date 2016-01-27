library(raster)
library(rasterVis)
library(ncdf4)

## READ THE NC

data <- stack("data.nc")

## valor máximo

maximos <- function(x){
    MM <- maxValue(max(x)) ## Encuentra el valor númerico máximo de la capa máxima del rasterStack
    xyMM <- xyFromCell(max(x), MM) ## Posición lon lat del máximo
    maxL <- which.max(x) ## calcula para cada celda la capa en la que se da el valor máximo
    Ly <- extract(maxL, xyMM) ## Para la posición de nuestro máximo extraemos el valor de la capa que nos dice cuándo s eda el máximo.
    sol <- c(MM, xyMM, Ly) ## Es el vector solución que contiene toda la información antes calculada.
    return(sol)
}

## Lo mismo para los mínimos

minimos <- function(x){
    mm <- minValue(min(x))
    xymm <- xyFromCell(min(x), mm)
    minL <- which.min(x)
    Ly <- extract(minL, xymm)
    sol <- c(mm, xymm, Ly)
    return(sol)
}

## Función para máximos y mínimos. Se selecciona lo que quieres calcular, bien máximos, bien mínimos en los argumentos de la función.

maxmin <- function(x, y) {
    if (y == 'max' ){
        MM <- maxValue(max(x))
        xyMM <- xyFromCell(max(x), MM)
        maxL <- which.max(x)
        Ly <- extract(maxL, xyMM)
        sol <- c(MM, xyMM, Ly)
        return(sol)
    }
    else {
        mm <- minValue(min(x))
        xymm <- xyFromCell(min(x), mm)
        minL <- which.min(x)
        Ly <- extract(minL, xymm)
        sol <- c(mm, xymm, Ly)
        return(sol)
    }
}
        
        
maximos(data)
minimos(data)
maxmin(data, 'max')
maxmin(data, 'min')





