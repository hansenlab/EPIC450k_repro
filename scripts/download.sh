# Esteller data: (GSE36369)
# Currently, the IDAT files are not available on GEO.
# We contacted the authors directly to obtain the IDAT files; 
# waiting for approval to release data or for the authors to add the IDAT files to their current GEO deposition

# In R: To download the package minfiDataEPIC containing the Illumina data:
source("https://bioconductor.org/biocLite.R")
biocLite("minfiDataEPIC")

# In R: To download the ENCODE 450k IDAT files from GEO:
library(GEOquery)
setwd("../data_raw/")
getGEOSuppFiles("GSE40699")