library(minfi)
library(minfiDataEPIC)
library(matrixStats)
library(pbapply)
source("../scripts/utils.R")
load("../data_processed/betas_only_epic.rda")


ann <- getAnnotation(RGsetEPIC)
probes <- getProbeType(RGsetEPIC, withColor=TRUE)
names(probes) <- rownames(ann)
probeNames <- rownames(ann)


# Reordeing:
betas <- pblapply(betas, function(x){
	x[probeNames, ]
})


vars <- do.call(cbind,pblapply(betas, rowVars))
rownames(vars) <- probeNames
types <- c("II", "IRed", "IGrn")
vars_list <- pblapply(types, function(x){
	indices <- probes==x
	vars[indices,]
})
names(vars_list) <- types

#Reordering:
names <- c("raw", "illumina", "swan", "quantile", "noob", "funnorm")
vars_list <- lapply(vars_list, function(x) x[, names])

# Seeting colors:
col_raw <- "grey70"
col_swan <- "red"
col_noob <- "orange"
col_illumina <- darken("olivedrab4",1.3)
col_funnorm <- "deeppink3"
col_quantile <- "deepskyblue3"


##########################################
##########################################
#############   Figure 1  ################
##########################################
##########################################
cols <- c(col_raw, col_illumina, col_swan, col_quantile, col_noob, col_funnorm)
labels <- c("Raw", "Illumina", "SWAN", "Quantile", "ssNoob", "Funnorm")
ylab="Probe variances"

pdf("../figures/probes_variances.pdf", width=9, height=6)
par(mfrow=c(1,2), bty="n")
#boxplot(vars_list[[2]], outline=FALSE)
boxplot(rbind(vars_list[[2]],vars_list[[3]]),
	outline=FALSE,
	main="Type I",
    ylab=ylab, 
    col=cols,
    xaxt="n",
    yaxt="n"
)
axis(side=2, at=c(0,0.0004, 0.0008, 0.0012))
text(1:6+0.2, 
	par("usr")[3]+0.5*par("usr")[3], 
	labels=labels,
	srt=45, 
	pos=2, 
	xpd=TRUE
)

boxplot(vars_list[[1]],
	outline=FALSE, 
	main="Type II", 
	ylab="", 
	col=cols,
	xaxt="n", 
	yaxt="n", 
	ylim=c(0, 0.0012)
)
axis(side=2, at=c(0,0.0004, 0.0008, 0.0012))
text(1:6+0.2, par("usr")[3]+0.5*par("usr")[3], 
	labels=labels,
	srt=45, 
	pos=2, 
	xpd=TRUE
)
dev.off()
##########################################
##########################################
##########################################
##########################################
##########################################




