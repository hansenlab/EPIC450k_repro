source("../scripts/utils.R")

extdir <- "../data_processed/"
load(file.path(extdir, "rgset_epic.rda"))
load(file.path(extdir, "rgset_450k.rda"))

epic_object_raw  <- preprocessRaw(rgset_epic)
epic_object_noob <- preprocessNoob(rgset_epic)
epic_object_illumina <- preprocessIllumina(rgset_epic)
epic_object_quantile <- preprocessQuantile(rgset_epic)
epic_object_swan     <- preprocessSWAN(rgset_epic)
epic_object_funnorm  <- preprocessFunnorm(rgset_epic)

m450k_object_raw   <- preprocessRaw(rgset_450k)
m450k_object_noob  <- preprocessNoob(rgset_450k)
m450k_object_illumina <- preprocessIllumina(rgset_450k)
m450k_object_quantile <- preprocessQuantile(rgset_450k)
m450k_object_swan     <- preprocessSWAN(rgset_450k)
m450k_object_funnorm  <- preprocessFunnorm(rgset_450k)

objects_epic <- list(raw=epic_object_raw, 
	illumina=epic_object_illumina, 
	swan=epic_object_swan, 
	quantile=epic_object_quantile, 
	noob=epic_object_noob, 
	funnorm=epic_object_funnorm
)

objects_450k <- list(raw=m450k_object_raw, 
	illumina=m450k_object_illumina, 
	swan=m450k_object_swan, 
	quantile=m450k_object_quantile, 
	noob=m450k_object_noob, 
	funnorm=m450k_object_funnorm
)

objects_combined <- pblapply(1:6, function(i){
	combineArrays(objects_450k[[i]], objects_epic[[i]])
})

betas <- pblapply(objects_combined, getBeta)
betas <- pblapply(betas, impute.matrix)
names(betas) <- names(objects_epic)
save(betas, file=file.path(extdir, "betas_combined_separate_normalization.rda"))