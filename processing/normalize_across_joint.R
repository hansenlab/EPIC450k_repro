source("../scripts/utils.R")

extdir <- "../extdata/epic_work"
load(file.path(extdir, "rgset_epic.rda"))
load(file.path(extdir, "rgset_450k.rda"))
rgset_combined <- combineArrayTypes(rgset_450k, rgset_epic)

object_raw  <- preprocessRaw(rgset_combined)
object_noob <- preprocessNoob(rgset_combined)
object_illumina <- preprocessIllumina(rgset_combined)
object_quantile <- preprocessQuantile(rgset_combined)
object_swan <- preprocessSWAN(rgset_combined)
object_funnorm <- preprocessFunnorm(rgset_combined)

objects <- list(raw=object_raw, 
	illumina=object_illumina, 
	swan=object_swan, 
	quantile=object_quantile, 
	noob=object_noob, 
	funnorm=object_funnorm
)

betas <- pblapply(objects, getBeta)
betas <- pblapply(betas, impute.matrix)
save(betas, file=file.path(extdir, "betas_combined.rda"))