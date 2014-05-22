####
library(sva)
library(limma)
library(scales)
library(grid)
library(RColorBrewer)

# this function regresses SVs out of genomic data
cleaningP = function(y, mod, svaobj,  P=ncol(mod)) {
        X=cbind(mod,svaobj$sv)
        Hat=solve(t(X)%*%X)%*%t(X)
        beta=(Hat%*%t(y))
        cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
        return(cleany)
}

# make a new directory to save plots
try(system("mkdir PCA_plots"))

# load data
load("normalizedES_GEO.rda")

p = theDataNorm$E
pd = theDataNorm$targets
map = theDataNorm$genes

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
