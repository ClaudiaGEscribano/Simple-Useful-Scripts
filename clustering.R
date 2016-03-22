## This script will do a Principal Component Analysis over a data set. After that kmeans algorithm will be implemented over the result.


## 1. Firts function: do the PCA for the data we provide.

pca2kmeans <- function(x){
    data <- as.data.frame(x)
    datapca <-princomp(data) ## PCA descomopsition
    ## I SELECT THE NUMBRE OF MODES THAT I CONSIDER. It is necessary to include in the function that this selection ewasne depending on a percentage.some conditions. This number (56) would be selected depending on some criteria.

    datapca2 <- data.frame(datapca$x[,1:56])
    return(datapca2)
}


datakm <- pca2kmeans(x)

## The data frame resulting from the function is the one we are going to use for the kmeans.

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

resultado <- kmeansexp(datakm,n,k) ## you can especify n=number of running times and k=number of clusters

## 3. 3rd function to analyse the correct number of cluters.

CH <- function(x, k){
    CH_index <- c()
    for (i in 2:k) CH_index[i] <- ((nrow(x)-i)/(i-1))*(x[[i]]$tot.withinss)/(x[[i]]$betweenss)
return(CH_index)
                               }

Indice <- lapply(resultado, FUN='CH')
Indice <- do.call(cbind, Indice) 
