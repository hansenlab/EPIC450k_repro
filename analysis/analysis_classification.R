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
		abs(x[,-epic_indices] - y)^2
	})
	# Creation of the rankings:
	medians <- pblapply(diffs, function(x){
		colMedians(x, na.rm=TRUE)
	})
	os <- pblapply(medians, order)
	list(os=os, medians=medians)
}


extdir <- "../data_processed/"
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
mains <- c("Raw", "Illumina", "SWAN", "Quantile", "ssNoob", "Funnorm")


ylim <- c(0,0.016)

pdf("../figures/classification_combined.pdf", width=5, height=7)
par(mfrow=c(6,1), mar=c(1,4,2,1))
for (i in 1:6){
	barplot(sort(os_combined$medians[[i]]), 
		col=col[os_combined$os[[i]]], 
		border=col[os_combined$os[[i]]], 
		ylim=ylim,
		main=mains[i], 
		yaxt="n",
		ylab="Median distance"
	)	
	axis(side=2, at=c(0, 0.008, 0.016), labels=c(0, 0.008, 0.016))
}
dev.off()

pdf("../figures/classification_separate.pdf", width=5, height=7)
par(mfrow=c(6,1), mar=c(1,4,2,1))
for (i in 1:6){
	barplot(sort(os_separate$medians[[i]]), 
		col=col[os_separate$os[[i]]], 
		border=col[os_separate$os[[i]]], 
		ylim=ylim,
		main=mains[i],
		yaxt="n",
		ylab="Median distance"
	)	
	axis(side=2, at=c(0, 0.008, 0.016), labels=c(0, 0.008, 0.016))
}
dev.off()



# Median distance for Gm12878
gm.index <- which(colnames(betas[[1]])=="Gm12878")
medians_combined <- unlist(lapply(betas_combined, function(beta){
	epic.mean <- rowMeans(beta[, 340:342])
	epic.gm   <- beta[, gm.index]
	diff <- (epic.mean-epic.gm)^2
	median(diff)
}))
medians_separate <- unlist(lapply(betas_separate, function(beta){
	epic.mean <- rowMeans(beta[, 340:342])
	epic.gm   <- beta[, gm.index]
	diff <- (epic.mean-epic.gm)^2
	median(diff)
}))



############################################################
####### # Dotplot
col_raw <- "black"
col_swan <- "red"
col_noob <- "orange"
col_illumina <- "olivedrab4"
col_quantile <- "deepskyblue3"
col_funnorm <- "deeppink3"
colors <- c(col_raw, col_illumina, col_swan, col_quantile, col_noob, col_funnorm)
names  <- c("Raw", "Illumina", "SWAN", "Quantile", "ssNoob", "Funnorm")
pdf("../figures/dotplot.pdf", width=4.5, height=5)
plot(medians_separate, 
	pch=20, cex=0.1, 
	bty="n", xaxt="n", 
	col="white", 
	ylim=c(0, 0.006),
	xlim=c(0, 6.5),
	xlab="",
	ylab="Median distance to ENCODE GM12878"
)
for (i in 1:6){
	abline(v=i, lty=3)
}
points(medians_combined, col=colors, cex=3, pch=20)
points(medians_separate, col=colors, cex=3)

text(1:6+0.35, 
	par("usr")[3]+1*par("usr")[3], 
	labels=names,
	srt=45, 
	pos=2, 
	xpd=TRUE
)
#legend(x=6, y=0.005, c("separate", "combined"), pch=c(1,20), cex=2, bty="n", title="Normalization")
dev.off()





############################################################
####### ROC CURVES
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