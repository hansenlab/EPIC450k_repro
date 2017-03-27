source("../scripts/utils.R")
extdir <- "../data_processed/"
load(file.path(extdir, "rgset_450k.rda"))
library(FlowSorted.Blood.450k)
blood <- rgset_450k[, 257:276]

blood.subset <- convertArray(blood, "IlluminaHumanMethylationEPIC")
blood.subset <- convertArray(blood.subset, "IlluminaHumanMethylation450k")

counts.full   <- estimateCellCounts(blood, returnAll=TRUE)
counts.subset <- estimateCellCounts(blood.subset, returnAll=TRUE)

# What is we convert to an EPIC array:
cors <- unlist(lapply(1:ncol(counts.full$counts), function(i){
	cor(counts.full$counts[,i], counts.subset$counts[,i])
}))
diff <- mean(abs(counts.full$counts - counts.subset$counts))