source("../scripts/utils.R")

extdir <- "../data_processed/"
load(file.path(extdir, "rgset_epic.rda"))


object_raw  <- preprocessRaw(rgset_epic)
object_noob <- preprocessNoob(rgset_epic)
object_swan <- preprocessSWAN(rgset_epic)
object_illumina <- preprocessIllumina(rgset_epic)
object_quantile <- preprocessQuantile(rgset_epic)
object_funnorm  <- preprocessFunnorm(rgset_epic)

objects <- list(raw=object_raw, 
	illumina= object_illumina,
	swan  = object_swan, 
	noob  = object_noob,
	quantile = object_quantile, 
	funnorm=object_funnorm
)

betas <- pblapply(objects, getBeta)
betas <- pblapply(betas, impute.matrix)
save(betas, file=file.path(extdir, "betas_only_epic.rda"))