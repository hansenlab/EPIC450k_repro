getOs <- function(betas){
	# Creating of a mean EPIC sample:
	means_epic <- pblapply(betas, function(x){
		epic_samples <- c("200144450018_R04C01","200144450019_R07C01","200144450021_R05C01")
		rowMeans(x[, epic_samples])
	})
	# Creation of the differences:
	diffs <- dlapply(betas, means_epic, function(x,y){
		epic_samples <- c("200144450018_R04C01","200144450019_R07C01","200144450021_R05C01")
		epic_indices <- match(epic_samples, colnames(x))
		abs(x[,-epic_indices] - y)
	})
	# Creation of the rankings:
	medians <- pblapply(diffs, function(x){
		colMedians(x, na.rm=TRUE)
	})
	os <- pblapply(medians, order)
	list(os=os, medians=medians)
}


extdir <- "../extdata/epic_work"
load(file.path(extdir, "betas_combined_separate_normalization.rda"))
betas_separate <- betas
load(file.path(extdir, "betas_combined.rda"))
betas_combined <- betas

# Measure one: variances statified by probe type: 
ann <- getAnnotation(RGsetEPIC)
probes <- getProbeType(RGsetEPIC, withColor=TRUE)
names(probes) <- rownames(ann)

# Getting ordered medians:
os_separate <- getOs(betas_separate)
os_combined <- getOs(betas_combined)



############################################################
####### Barplots

# Pheno and color stuff:
pheno <- rep("LCL (Esteller)", ncol(betas[[1]])-3)
pheno[277:339] <- "ENCODE"
pheno[257:276] <- "PBMC"
pheno[match(c("Gm19239", 
	"Gm06990", 
	"Gm12891", 
	"Gm12878", 
	"Gm12892"), colnames(betas[[1]]))] <- "LCL (ENCODE)"
table(pheno)

col1 <- col11 <- "grey"
col2 <- "deepskyblue3"
col3 <- "firebrick"
col <- rep(NA, length(pheno))
col[pheno=="LCL (Esteller)"] <- col1
col[pheno=="LCL (ENCODE)"]   <- col11
col[pheno=="PBMC"] <- col2
col[pheno=="ENCODE"] <- col3
mains <- c("Raw", "Illumina", "SWAN", "Quantile", "noob", "Funnorm")


pdf("../figures/classification_combined.pdf", width=5, height=7)
par(mfrow=c(6,1), mar=c(1,4,2,1))
for (i in 1:6){
	barplot(sort(os_combined$medians[[i]]), 
		col=col[os_combined$os[[i]]], 
		border=col[os_combined$os[[i]]], 
		ylim=c(0,0.12), 
		main=mains[i], 
		yaxt="n",
		ylab="Median distance"
	)	
	axis(side=2, at=c(0, 0.06, 0.12), labels=c(0, 0.06, 0.12))
}
dev.off()

pdf("../figures/classification_separate.pdf", width=5, height=7)
par(mfrow=c(6,1), mar=c(1,4,2,1))
for (i in 1:6){
	barplot(sort(os_separate$medians[[i]]), 
		col=col[os_separate$os[[i]]], 
		border=col[os_separate$os[[i]]], 
		ylim=c(0,0.12),
		main=mains[i],
		yaxt="n",
		ylab="Median distance"
	)	
	axis(side=2, at=c(0, 0.06, 0.12), labels=c(0, 0.06, 0.12))
}
dev.off()




############################################################
####### ROC CURVES
col_raw <- "black"
col_swan <- "red"
col_noob <- "orange"
col_illumina <- "olivedrab4"
col_quantile <- "deepskyblue3"
col_funnorm <- "deeppink3"
colors <- c(col_raw, col_illumina, col_swan, col_quantile, col_noob, col_funnorm)
names  <- c("Raw", "Illumina", "SWAN", "Quantile", "noob", "Funnorm")

# Binary outcome for ROC curves:
y <- rep(0, length(pheno))
y[grepl("LCL", pheno)] <- 1


pdf("../figures/roc_combined.pdf", width=5, height=5)
par(bty="n")
lwd=2
plot(roc(y, os_combined$medians[[1]]), col=col_raw, lwd=lwd)
for (i in 2:6){
	lines(roc(y, os_combined$medians[[i]]), col=colors[i], lwd=lwd)
}
legend("bottomright", names, col=colors, lwd=lwd, bty="n")
dev.off()


pdf("../figures/roc_separate.pdf", width=5, height=5)
par(bty="n")
lwd=2
plot(roc(y, os_separate$medians[[1]]), col=col_raw, lwd=lwd)
for (i in 2:6){
	lines(roc(y, os_separate$medians[[i]]), col=colors[i], lwd=lwd)
}
legend("bottomright", names, col=colors, lwd=lwd, bty="n")
dev.off()