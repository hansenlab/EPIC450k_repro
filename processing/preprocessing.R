# This script builds the RGChannelSets from the raw idat files.
# This script has to be run after the data are downloaded with the script download.sh
source("../scripts/utils.R")


#### PROCESSING ENCODE:
files <- unique(gsub("_Grn.idat|_Red.idat","",list.files("../data_raw/GSE40699/", pattern=".idat")))
rgset_encode <- read.450k(files)
colnames(rgset_encode) <- gsub("GSM[0-9]*_hg19_wgEncodeHaibMethyl450|SitesRep1","",colnames(rgset_encode))
sampleNames(phenoData(rgset_encode)) <- sampleNames(assayData(rgset_encode))
rgset_encode <- updateObject(rgset_encode)
save(rgset_encode, file="../data_processed/rgset_encode.rda")


#### PROCESSING ESTELLER:
# At the moment, the IDAT files are not publicly available
sheet <- read.csv("/dcl01/hansen/data/esteller_lcl/HVP/HVP.csv", skip=7, head=TRUE)
files <- list.files("/dcl01/hansen/data/esteller_lcl/HVP", 
	recursive=TRUE, pattern="idat", full.names=TRUE)
files <- files[!duplicated(gsub("Grn|Red", "", files))]
rgset <- read.450k(files)

sheet$filename <- paste0(sheet$Sentrix_ID,"_", sheet$Sentrix_Position)
sheet <- sheet[match(colnames(rgset), sheet$filename),]
col <- as.character(sheet$Sample_Name)
col <- substr(col,1,2)
col <- as.numeric(as.factor(col))
sheet$group <- col
pData(rgset) <- sheet
rgset_esteller <- rgset
sampleNames(phenoData(rgset_esteller)) <- sampleNames(assayData(rgset_esteller))
rgset_esteller <- updateObject(rgset_esteller)
save(rgset_esteller, file="../data_processed/rgset_esteller.rda")


# Building the 450k RGChannelSet
rgset_450k <- combine(rgset_esteller, rgset_encode)
pData(rgset_450k) <- pData(rgset_450k)[,1, drop=FALSE]
pData_450k <- pData(rgset_450k)
pData_450k$Sample_Name <- as.character(pData_450k$Sample_Name) 
pData(rgset_450k) <- pData_450k
save(rgset_450k, file="../data_processed/rgset_450k.rda")

# Building the EPIC RGChannelSet
rgset_epic <- RGsetEPIC
pData(rgset_epic) <- pData(rgset_epic)[,1, drop=FALSE]
pData_epic <- pData(rgset_epic)
pData_epic$Sample_Name <- as.character(pData_epic$Sample_Name) 
pData(rgset_epic) <- pData_epic
save(rgset_epic, file="../data_processed/rgset_epic.rda")