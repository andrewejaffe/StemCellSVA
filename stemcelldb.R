#########
# Preprocess StemCellDB in GEO
# Andrew Jaffe
# Updated 4/8/13
#########################

##### INSTRUCTIONS ######
# 1. PLACE "G4112F_annotation.txt" IN YOUR WORKING DIRECTORY
# 2. SOURCE FUNCTIONS AT THE TOP OF THIS SCRIPT (SEARCH TERM: 'FXN')
# 3. PREPROCESSING (SEARCH TERM: 'PREPROC')
# 4. MAKE PCA FIGURES (SEARCH TERM: 'FIG1')
# 5. MAKE D.E. FIGURES (SEARCH TERM: 'FIG2')
# 6. MAKE SEX AND GSTT1 FIGURES (SEARCH TERM: 'FIG3')

###### LOAD LIBRARIES ############
# if these aren't installed, here's the BioC install:
# source("http://bioconductor.org/biocLite.R")

## BIOC
library(limma)
library(GEOquery)
library(sva)

## CRAN
library(scales)
library(RColorBrewer)


########################################
#### FXN: SOME FUNCTIONS
# splits character strings, wrapper for strsplit
ss=function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), function(y) y[slot])

# this function regresses SVs out of genomic data
cleaningP = function(y, mod, svaobj,  P=ncol(mod)) {
        X=cbind(mod,svaobj$sv)
        Hat=solve(t(X)%*%X)%*%t(X)
        beta=(Hat%*%t(y))
        cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
        return(cleany)
}

# splits in factors
splitit=function(x) split(seq(along=x),x)


#####
# plotting code, expBySample()
# pp is a vector of expression values, or 1 row, of the expression matrix
# sampleIds is a vector of sample IDs
# conditions: colors the samples by treatment or condition, a vector
expBySample = function(pp, sampleIds, geneName, conditions=NULL,
	pal="Dark2",legend=TRUE,cex=0.85,...) {
	
	require(RColorBrewer)
	
	sampleIds = factor(sampleIds)
	N=length(levels(sampleIds))
	if(!is.null(conditions)) {
		conditions = factor(conditions)
		cols = as.numeric(conditions)
		Ncol = max(cols)
	} else {
		cols = rep(1,length(pp))
		Ncol = 1
	}
	
	if(pal %in% rownames(brewer.pal.info)) {
		palette(brewer.pal(max(Ncol,3), pal))
	} else palette(pal)

	xx = jitter(as.numeric(sampleIds),amount=0.25)
	plot(pp~xx,	pch = 21, bg=cols, cex=cex,
		xaxt = "n", ylab = "Expression", xlab="", 
		main = geneName,...)
	abline(v=seq(0.5,N+0.5, by=1),lty=2, lwd = 0.8)	
	axis(1, at = 1:N, labels = levels(sampleIds), las = 2,...)
	if(legend) legend("topright", levels(conditions),
				col = 1:Ncol,pch=15, nc =Ncol,bty="n",cex=0.9)
}

###############
## PREPROC: PREPROCESSING 
tmp = getGEO("GSE32923")
pd = pData(tmp[[1]])
pd = pd[,c(1,2, grep("characteristics",names(pd)))]

# CLEAN UP PHENOTYPE DATA
for(i in 3:ncol(pd)) {
	names(pd)[i] = ss(as.character(pd[,i]), ": ", 1)[1]
	pd[,i] = ss(as.character(pd[,i]), ": ", 2)
}
names(pd) = c("Title", "GEO_Accession", "CellType",	
	"Treatment", "SampleID", "Passage", "Sex")
	
## only keep ESC and iPSC, drop last column
pd = pd[pd$CellType %in% c("ESC","iPSC"),-ncol(pd)]

## raw data
tmp = getGEOSuppFiles("GSE32923")
system("tar xvf GSE32923/GSE32923_RAW.tar -C GSE32923/")
system("gunzip GSE32923/*.gz GSE32923/")

### match to array data
id = ss(dir("GSE32923", pattern = "GSM"), "_",1)
fns = dir("GSE32923", pattern = "GSM",full.names=TRUE)
pd$FileName = fns[match(pd$GEO_Accession, id)]
pd$Passage = as.numeric(pd$Passage)

## read in data
theData = read.maimages(pd$FileName, source="agilent",
	green.only=TRUE, columns = list(G = "gMedianSignal", Gb = "gBGMedianSignal"),
	annotation = c("accessions","chr_coord", "ControlType",
	"ProbeName", "GeneName", "SystematicName", "Description","Sequence"))

# get processing date
pd$ArrayDate = sapply(pd$FileName, function(x) {
	tmp = read.delim(file=x,skip=1,nrows=1,header=T,as.is=T)
	out = strsplit(tmp[1,"FeatureExtractor_ExtractionTime"], " ")[[1]][1]
	return(out)
})
theData$targets = pd

######## PREPROC  #########

# drop control probes
keepIndex = which(theData$genes$ControlType == 0)
theDataClean = theData[keepIndex,]

# for now, default background correct
theDataBack = backgroundCorrect(theDataClean, offset = 50)

# lets normalize all cell types together for now
theDataNorm = normalizeBetweenArrays(theDataClean, method="quantile")

#################
# more annotation
map = theDataNorm$genes

### UPDATE ANNOTATION, THIS TXT FILE IS ON GITHUB
tmp = read.delim("G4112F_annotation.txt",	header = T, as.is=T)
mIndex = match(map$ProbeName, tmp$ProbeID)
map = tmp[mIndex,]
map$probeLoc = theDataNorm$genes$chr_coord
map$ProbeSequence = theDataNorm$genes$Sequence
rownames(theDataNorm$E) = map$ProbeID

map$geneChr = sapply(strsplit(map$GenomicCoordinates,":"), function(x) x[1])
pos = sapply(strsplit(map$GenomicCoordinates,":"), function(x) x[2])
map$geneStart = sapply(strsplit(pos,"-"), function(x) x[1])
map$geneEnd = sapply(strsplit(pos,"-"), function(x) x[2])

map$geneStrand = ifelse(map$geneStart < map$geneEnd, "+","-")

theDataNorm$genes = map

# save
save(theDataNorm, file = "normalizedES_GEO.rda",
	compress=TRUE)
	
#########################
## FIG1: FIGURE 1 - PCA
########################

# make a new directory to save plots
try(system("mkdir PCA_plots"))

p = theDataNorm$E
pd = theDataNorm$targets
map = theDataNorm$genes

#update sex
sexIndex = grep("DDX3Y",map$GeneSymbol)
sexCheck = colMeans(p[sexIndex,]) 
pd$Sex = ifelse(sexCheck > 9, "Male","Female")

##### plot raw data
pd$dates = factor(as.Date(pd$ArrayDate, format="%d-%b-%Y"))

# DO PCA
pca = prcomp(t(p))
pcaVars=signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100
signed = ifelse(max(pca$x[,2] > 70), 1, -1) # same sign across plots

#####
# Figure 1A: pca on ES, colored by treatment
pdf("PCA_plots/figure_1A_ES_PCA_noSVA.pdf")
palette(brewer.pal(8,"Dark2"))
 
plot(pca$x[,1], signed*pca$x[,2], 
	pch = 19,col = as.numeric(factor(pd$Treatment)),
	main = "Figure 1A - No SVA", cex=1.1,
	ylim = c(-80,100),
	cex.axis=1.5,cex.lab=1.5,
	xlab = paste("PC1:",pcaVars[1],"% of Variance Explained"),
	ylab = paste("PC2:",pcaVars[2],"% of Variance Explained"))
legend("bottom", c("FSB","KSR","UNDIFF"), col = 1:3, 
	pt.cex = 2, nc = 3, cex=1.1, pch = 15)
dev.off()

##########
# Figure 1B: No SVA, colored by Batch
n = length(levels(pd$dates))

batchNum = c(letters,LETTERS)[as.numeric(pd$dates)]
colPal = colorRampPalette(brewer.pal(9,"YlOrRd")[2:9])

pdf("PCA_plots/figure_1B_batch_noSVA.pdf")
palette(colPal(n))
plot(pca$x[,1], signed*pca$x[,2], 
	type="n",	main = "", cex=1.1,
	ylim = c(-80,100),	cex.axis=1.5,cex.lab=1.5,
	xlab = paste("PC1:",pcaVars[1],"% of Variance Explained"),
	ylab = paste("PC2:",pcaVars[2],"% of Variance Explained"))
text(pca$x[,1], signed*pca$x[,2], batchNum,
	col="black",font=2, cex=1.05, adj = c(0.5,0.5))
text(pca$x[,1], signed*pca$x[,2], batchNum,
	col=as.numeric(pd$dates), cex=0.96,adj = c(0.5,0.5))
dev.off()

# inset
pdf("PCA_plots/figure_1B_inset_pc2_vs_batch.pdf",	w=4)
palette(colPal(n))
boxplot(signed*pca$x[,2]~pd$dates,	horiz = T,
	ylab = "",	xaxt="n", ylim = c(-80,100),
	col=1:n)
	
lets = c(letters,LETTERS)[1:n]
text(seq(1,n,by=2), y = 95, lets[seq(1,n,by=2)])
text(seq(2,n,by=2), y = 100, lets[seq(2,n,by=2)])
text(1:n, -75, levels(pd$dates), cex=0.6, srt="270",font=2)
dev.off()

########################
# SVA

mod = model.matrix(~Treatment, data =pd)

# sva with default # of SVs
svaobj = sva(p, mod, n.sv=27) # auto is 27

# regress out SVs
cleanp = cleaningP(p,mod,svaobj)

# run PCA
pca2 = prcomp(t(cleanp))
pcaVars2=signif(((pca2$sdev)^2)/(sum((pca2$sdev)^2)),3)*100
signed = ifelse(max(pca2$x[,2] > 70), 1, -1) # make same sign

#####
# Figure 1C: pca on ES, colored by treatment
pdf("PCA_plots/figure_1C_ES_PCA_SVA.pdf")
palette(brewer.pal(8,"Dark2"))
plot(pca2$x[,1], signed*pca2$x[,2], 
	pch = 19,col = as.numeric(factor(pd$Treatment)),
	main = "Figure 1C - SVA", cex=1.1,
	ylim = c(-35,65),
	cex.axis=1.5,cex.lab=1.5,
	xlab = paste("PC1:",pcaVars2[1],"% of Variance Explained"),
	ylab = paste("PC2:",pcaVars2[2],"% of Variance Explained"))
dev.off()

##########
# Figure 1D: SVA, colored by Batch

pdf("PCA_plots/figure_1D_batch_SVA.pdf")
palette(colPal(n))
plot(pca2$x[,1], signed*pca2$x[,2], 
	type="n",	main = "", cex=1.1,
	ylim = c(-35,65),cex.axis=1.5,cex.lab=1.5,
	xlab = paste("PC1:",pcaVars2[1],"% of Variance Explained"),
	ylab = paste("PC2:",pcaVars2[2],"% of Variance Explained"))
text(pca2$x[,1], signed*pca2$x[,2], batchNum,
	col="black",font=2, cex=1.05, adj = c(0.5,0.5))
text(pca2$x[,1], signed*pca2$x[,2], batchNum,
	col=as.numeric(pd$dates), cex=0.96,adj = c(0.5,0.5))
dev.off()

# inset
pdf("PCA_plots/figure_1D_inset_pc2_vs_batch_SVA.pdf",	w=4)
palette(colPal(n))
boxplot(signed*pca2$x[,2]~pd$dates, horiz = T,
	ylab = "",	xaxt="n",ylim = c(-35,65),
	col=1:n)
	
lets = c(letters,LETTERS)[1:n]
text(seq(1,n,by=2), y = 65, lets[seq(1,n,by=2)])
text(seq(2,n,by=2), y = 62, lets[seq(2,n,by=2)])
text(1:n, -33, levels(pd$dates), cex=0.6, srt="270", font=2)
dev.off()

#### 
# FIG2: DIFFERENTIAL EXPRESSION

# fit linear model without SVA
fit0 = lmFit(p, model.matrix(~pd$Treatment))
eb0 = ebayes(fit0) # empirical bayes for t-stats

# with SVA
fit1 = lmFit(p, cbind(model.matrix(~pd$Treatment), svaobj$sv))
eb1 = ebayes(fit1)

### make a new folder for plots
try(system("mkdir treat_fig"))

####
# plot treatment effects
oIndex = grep("OLFML1", map$GeneSymbol)

pdf("treat_fig/treat_es_noSVA_fig2a.pdf")
expBySample(p[oIndex,], pd$SampleID, geneName = "",cex=1,
	conditions = pd$Treatment,legend = FALSE, cex.axis=1.3, cex.lab=1.5)
dev.off()

pdf("treat_fig/treat_es_SVA_fig2b.pdf")
expBySample(cleanp[oIndex,], pd$SampleID, geneName = "",cex=1,
	conditions = pd$Treatment,legend = FALSE, cex.axis=1.3, cex.lab=1.5)
legend("bottomright", levels(factor(pd$Treatment)), bty="n",
	col = 1:3, pch = 15, pt.cex = 3,cex=1.6, nc = 3)
dev.off()

#############
# check pax6
pIndex = grep("PAX6", map$GeneSymbol)
eb0$p[pIndex,2:3]
eb1$p[pIndex,2:3]

pdf("treat_fig/treat_es_noSVA_pax6.pdf")
expBySample(p[pIndex,], pd$SampleID, geneName = "",cex=1,
	conditions = pd$Treatment,legend = FALSE, cex.axis=1.3, cex.lab=1.5)
dev.off()

pdf("treat_fig/treat_es_SVA_pax6.pdf")
expBySample(cleanp[pIndex,], pd$SampleID, geneName = "",cex=1,
	conditions = pd$Treatment,legend = FALSE, cex.axis=1.3, cex.lab=1.5)
dev.off()

#### SUPP FIGURE
# plot p-value
pdf("treat_fig/overall_pval.pdf")
plot(-log10(eb1$p[,2]), -log10(eb0$p[,2]), 
	xlab="p-value (SVA)",ylab="p-value (No SVA)",
	col = alpha("black",0.33),
	cex.axis=1.5)
points(-log10(eb1$p[,2])[oIndex], -log10(eb0$p[,2])[oIndex],
	col = "red",pch=19,cex=1.4)
points(-log10(eb1$p[,2])[pIndex], -log10(eb0$p[,2])[pIndex],
	col = "blue",pch=19,cex=1.4)
abline(0,1,col="blue")
dev.off()

# plot sigma
pdf("treat_fig/overall_sigma.pdf")
plot(fit1$sigma, fit0$sigma, 
	xlab="sigma (SVA)",ylab="sigma (No SVA)",
	col = alpha("black",0.33),
	cex.axis=1.5)
abline(0,1,col="blue")
points(fit1$sigma[oIndex], fit0$sigma[oIndex],
	col = "red",pch=19,cex=1.4)
points(fit1$sigma[pIndex], fit0$sigma[pIndex],
	col = "blue",pch=19,cex=1.4)
dev.off()

# plot betas
pdf("treat_fig/overall_betas.pdf")
plot(fit1$coef[,2], fit0$coef[,2], 
	xlab="beta (SVA)",ylab="beta (No SVA)",
	col = alpha("black",0.33),
	cex.axis=1.5)
abline(0,1,col="blue")
points(fit1$coef[oIndex,2], fit0$coef[oIndex,2],
	col = "red",pch=19,cex=1.4)
points(fit1$coef[pIndex,2], fit0$coef[pIndex,2],
	col = "blue",pch=19,cex=1.4)
dev.off()


##############################
# FIG3: GSTT1 PLOTS

try(system("mkdir sex_gstt1"))

gIndex = grep("GSTT1", map$GeneSymbol)[2]

ylim = c(6,12.5)
pal = c("black","black")

# take the best hit
pdf("sex_gstt1/gstt1_es_noSVA_fig2d.pdf")
expBySample(p[gIndex,], pd$SampleID, geneName = "",cex=1,ylim=ylim,
	conditions = NULL, pal=pal,legend = FALSE, cex.axis=1.3, cex.lab=1.5)
dev.off()

pdf("sex_gstt1/gstt1_es_wrongSVA_fig2e.pdf")
expBySample(cleanp[gIndex,], pd$SampleID, geneName = "",cex=1,
	conditions = NULL, pal=pal,legend = FALSE, cex.axis=1.3, cex.lab=1.5,
	ylim = ylim)
dev.off()

# determine CNV by group
gstt1exp = tapply(p[gIndex,], pd$SampleID, median)
plot(gstt1exp)
gstt1pred = cut(gstt1exp, c(5,7.5, 10.5,12), labels = F)-1
names(gstt1pred) = names(gstt1exp)

pd$gstt1 = rep(NA)
sIndexes = splitit(pd$SampleID)
for(i in seq(along=sIndexes)) {
	ind=sIndexes[[i]]
	pd$gstt1[ind] = rep(gstt1pred[names(sIndexes)[i]])
}

# SVA FOR GSTT1 MODEL
mod3 = model.matrix(~gstt1+Treatment, data =pd)
svaobj3 = sva(p,mod3, n.sv=26) # default
cleanp3 = cleaningP(p,mod3,svaobj3)

pdf("sex_gstt1/gstt1_es_rightSVA_fig2f.pdf")
expBySample(cleanp3[gIndex,], pd$SampleID, geneName = "",cex=1,
	conditions = NULL, pal=pal,legend = FALSE, cex.axis=1.3, cex.lab=1.5,
	ylim = ylim)
dev.off()

#########
# FIG3: SEX DIFFERENCES

# find sex difference w/o SVA
fit0a = lmFit(p, model.matrix(~pd$Sex + pd$Treatment))
eb0a = ebayes(fit0a)
sIndex = grep("RPS4Y1", map$GeneSymbol)


# take the best hit
pdf("sex_gstt1/sex_es_noSVA_fig2a.pdf")
expBySample(p[sIndex,], pd$SampleID, geneName = map$GeneSymbol[sIndex],
	conditions = pd$Sex, legend = FALSE, cex.axis=1.3, cex.lab=1.5,cex=1,
	pal="Set1")
dev.off()

pdf("sex_gstt1/sex_es_wrongSVA_fig2b.pdf")
expBySample(cleanp[sIndex,], pd$SampleID, geneName = map$GeneSymbol[sIndex],
	conditions = pd$Sex, legend = FALSE, cex.axis=1.3, cex.lab=1.5,
	ylim = range(p[sIndex,]), cex=1, pal="Set1")
legend("top", levels(factor(pd$Sex)), col = 1:2, pch = 15, 
	pt.cex=3, nc = 2, cex =1.75,bty="n")
dev.off()


#### adjust for sex in SVA
mod2 = model.matrix(~Sex+Treatment, data =pd)
svaobj2 = sva(p,mod2, n.sv=26)
cleanp2 = cleaningP(p,mod2,svaobj2)

pdf("sex_gstt1/sex_es_rightSVA_fig2c.pdf")
expBySample(cleanp2[sIndex,], pd$SampleID, geneName = map$GeneSymbol[sIndex],
	conditions = pd$Sex, legend = FALSE, cex.axis=1.3, cex.lab=1.5,
	ylim = range(p[sIndex,]), cex=1, pal="Set1")
dev.off()
