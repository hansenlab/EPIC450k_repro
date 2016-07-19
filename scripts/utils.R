library(pbapply)
library(matrixStats)
library(pROC)
library(minfi)
library(minfiDataEPIC)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

impute.matrix <- function (matrix) {
	missing <- which(is.na(matrix) | !is.finite(matrix), arr.ind=TRUE)
	if (length(missing)!=0){
		for (j in 1:nrow(missing)){
				mean <- mean(matrix[missing[j,1],][is.finite(matrix[missing[j,1],])], na.rm=TRUE)
				matrix[missing[j,1],missing[j,2]] <- mean
			}
	}
  	matrix
}

darken <- function(color, factor=1.4){
    col <- col2rgb(color)
    col <- col/factor
    col <- rgb(t(col), maxColorValue=255)
    col
}

dlapply <- function(X,Y, FUN){
	temp <- list()
	n <- length(X)
	for (i in 1:n){
		temp[[i]] <- FUN(X[[i]], Y[[i]])
	}
	temp
}