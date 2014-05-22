####
library(sva)
library(limma)
library(scales)
library(RColorBrewer)

# this function regresses SVs out of genomic data
cleaningP = function(y, mod, svaobj,  P=ncol(mod)) {
        X=cbind(mod,svaobj$sv)
        Hat=solve(t(X)%*%X)%*%t(X)
        beta=(Hat%*%t(y))
        cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
        return(cleany)
}
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

## load the data
load("normalizedES_GEO.rda")

p = theDataNorm$E
pd = theDataNorm$targets
map = theDataNorm$genes

####
pd$Treatment = factor(pd$Treatment, levels = c("UNDIFF", "FBS", "KSR"))

mod = model.matrix(~Treatment, data =pd)
svaobj = sva(p,mod, n.sv=27) # do SVA, 27 is default
cleanp = cleaningP(p,mod,svaobj) # cleaned data

######
# find treatment differences w and w/o SVA

# fit linear model without SVA
fit0 = lmFit(p, mod)
eb0 = ebayes(fit0) # empirical bayes for t-stats

# with SVA
fit1 = lmFit(p, cbind(mod, svaobj$sv))
eb1 = ebayes(fit1)

### make a new folder for plots
try(system("mkdir treat_fig"))

## other known genes
neuroGeneSymbols = c("TSPAN31", "BBS10", "TARBP1","C17orf69",
	"ABI2", "NRCAM", "SCG5", "SMAD7", "CORO2A", "SLC25A24")
mesenGeneSymbols = c("MBNL2", "ZHX1", "TOM1L2", "HGSNAT", "SSFA2", 
	"FKBP9", "PARVA",
	"LRRC8B", "NPM3", "SNURF", "RTN3", "EIF5A2", "CUTC", "SALL2")
neurIndex = match(neuroGeneSymbols, map$GeneSymbol)
mesenIndex = match(mesenGeneSymbols, map$GeneSymbol)

####
# plot treatment effects
oIndex = grep("OLFML1", map$GeneSymbol)
pd$TreatmentColor = factor(pd$Treatment)

pdf("treat_fig/treat_es_noSVA_fig2a_check.pdf")
expBySample(p[oIndex,], pd$SampleID, geneName = "",cex=1,
	conditions = as.character(pd$Treatment),
	legend = FALSE, cex.axis=1.3, cex.lab=1.5)
dev.off()

pdf("treat_fig/treat_es_SVA_fig2b.pdf")
expBySample(cleanp[oIndex,], pd$SampleID, geneName = "",cex=1,
	conditions = as.character(pd$Treatment),
	legend = FALSE, cex.axis=1.3, cex.lab=1.5)
legend("bottomright", levels(factor(as.character(pd$Treatment))), bty="n",
	col = 1:3, pch = 15, pt.cex = 3,cex=1.6, nc = 3)
dev.off()

plot(eb1$t[neurIndex,3], eb0$t[neurIndex,3],
	main = "neuro genes", xlab="SVA", ylab="No SVA",
	ylim = c(-20,20), xlim=c(-20,20))
abline(0,1)
plot(rank(eb1$p[,3])[neurIndex], rank(eb0$p[,3])[neurIndex])
abline(0,1)

sort(rank(eb1$p[,3])[neurIndex])

pdf("suppFigure_ranks_neuro.pdf")
plot(rank(eb1$p[,3]), rank(eb0$p[,3]),pch=21,bg="grey",
	xlab="Post-SVA Rank", ylab="Pre-SVA Rank",cex=0.7,
	cex.axis=1.5,cex.lab=1.5)
points(rank(eb1$p[,3])[neurIndex], cex=2,
	rank(eb0$p[,3])[neurIndex],pch=21,bg="red")
abline(0,1,col="blue",lty=2)	
dev.off()

pdf("suppFigure_ranks_mesen.pdf")
plot(rank(eb1$p[,2]), rank(eb0$p[,2]),pch=21,bg="grey",
	xlab="Post-SVA Rank", ylab="Pre-SVA Rank",cex=0.7)
points(rank(eb1$p[,2])[mesenIndex], cex=2,
	rank(eb0$p[,2])[mesenIndex],pch=21,bg="red")
abline(0,1,col="blue",lty=2)	

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

# plot p-value
pdf("treat_fig/overall_pval_new.pdf")
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
# gstt1

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

mod3 = model.matrix(~gstt1+Treatment, data =pd)
svaobj3 = sva(p,mod3, n.sv=26) # default
cleanp3 = cleaningP(p,mod3,svaobj3)

pdf("sex_gstt1/gstt1_es_rightSVA_fig2f.pdf")
expBySample(cleanp3[gIndex,], pd$SampleID, geneName = "",cex=1,
	conditions = NULL, pal=pal,legend = FALSE, cex.axis=1.3, cex.lab=1.5,
	ylim = ylim)
dev.off()

#########
# sex

# find sex difference w/o SVA
fit0a = lmFit(p, model.matrix(~pd$Sex + pd$Treatment))
eb0a = ebayes(fit0a)
sIndex = grep("RPS4Y1", map$GeneSymbol)


# take the best hit
pdf("sex_gstt1/sex_es_noSVA_fig2a.pdf")
expBySample(p[sIndex,], pd$SampleID, geneName = map$GeneSymbol[sIndex],
	conditions = pd$Sex, legend = FALSE, cex.axis=1.3, cex.lab=1.5,cex=1,
	ylim = c(4,17), pal="Set1")
dev.off()

pdf("sex_gstt1/sex_es_wrongSVA_fig2b.pdf")
expBySample(cleanp[sIndex,], pd$SampleID, geneName = map$GeneSymbol[sIndex],
	conditions = pd$Sex, legend = FALSE, cex.axis=1.3, cex.lab=1.5,
	ylim = c(4,17), cex=1, pal="Set1")
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
	ylim = c(4,17), cex=1, pal="Set1")
dev.off()
