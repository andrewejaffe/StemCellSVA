###
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

# load data
load("normalizedES_GEO.rda")

p = theDataNorm$E
pd = theDataNorm$targets
map = theDataNorm$genes

mod = model.matrix(~Treatment, data =pd)
batch = pd$SampleID

cleanp = ComBat(p, batch, mod)

# run PCA
pca2 = prcomp(t(cleanp))
pcaVars2=signif(((pca2$sdev)^2)/(sum((pca2$sdev)^2)),3)*100
signed = ifelse(max(pca2$x[,2] > 70), 1, -1) # make same sign

#####
# Figure 1C: pca on ES, colored by treatment
pdf("PCA_plots/supp_figure_ComBat_A.pdf")
palette(brewer.pal(8,"Dark2"))
plot(pca2$x[,1], signed*pca2$x[,2], 
	pch = 19,col = as.numeric(factor(pd$Treatment)),
	main = "Figure 1C - ComBat", cex=1.1,
	ylim = c(-35,65),
	cex.axis=1.5,cex.lab=1.5,
	xlab = paste("PC1:",pcaVars2[1],"% of Variance Explained"),
	ylab = paste("PC2:",pcaVars2[2],"% of Variance Explained"))
dev.off()

##########
# Figure 1D: SVA, colored by Batch
pd$dates = factor(as.Date(pd$ArrayDate, format="%d-%b-%Y"))
n = length(levels(pd$dates))
batchNum = c(letters,LETTERS)[as.numeric(pd$dates)]
colPal = colorRampPalette(brewer.pal(9,"YlOrRd")[2:9])

pdf("PCA_plots/supp_figure_ComBat_B.pdf")
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
pdf("PCA_plots/supp_figure_ComBat_C.pdf",	w=4)
palette(colPal(n))
boxplot(signed*pca2$x[,2]~pd$dates, horiz = T,
	ylab = "",	xaxt="n",ylim = c(-35,65),
	col=1:n)
	
lets = c(letters,LETTERS)[1:n]
text(seq(1,n,by=2), y = 65, lets[seq(1,n,by=2)])
text(seq(2,n,by=2), y = 62, lets[seq(2,n,by=2)])
text(1:n, -33, levels(pd$dates), cex=0.6, srt="270", font=2)
dev.off()

################################
### combat wrong model
mod2 = model.matrix(~1, data=pd)
cleanp2 = ComBat(p, batch, mod2)


# run PCA
pca3 = prcomp(t(cleanp2))
pcaVars3=signif(((pca3$sdev)^2)/(sum((pca3$sdev)^2)),3)*100
signed = ifelse(max(pca3$x[,2] > 70), 1, -1) # make same sign

#####
# Figure 1C: pca on ES, colored by treatment
pdf("PCA_plots/supp_figure_ComBat_D.pdf")
palette(brewer.pal(8,"Dark2"))
plot(pca3$x[,1], signed*pca3$x[,2], 
	pch = 19,col = as.numeric(factor(pd$Treatment)),
	main = "Figure 1D - ComBat", cex=1.1,
	ylim = c(-35,65),
	cex.axis=1.5,cex.lab=1.5,
	xlab = paste("PC1:",pcaVars3[1],"% of Variance Explained"),
	ylab = paste("PC2:",pcaVars3[2],"% of Variance Explained"))
dev.off()

pdf("PCA_plots/supp_figure_ComBat_E.pdf")
palette(colPal(n))
plot(pca3$x[,1], signed*pca3$x[,2], 
	type="n",	main = "", cex=1.1,
	ylim = c(-35,65),cex.axis=1.5,cex.lab=1.5,
	xlab = paste("PC1:",pcaVars3[1],"% of Variance Explained"),
	ylab = paste("PC2:",pcaVars3[2],"% of Variance Explained"))
text(pca3$x[,1], signed*pca3$x[,2], batchNum,
	col="black",font=2, cex=1.05, adj = c(0.5,0.5))
text(pca3$x[,1], signed*pca3$x[,2], batchNum,
	col=as.numeric(pd$dates), cex=0.96,adj = c(0.5,0.5))
dev.off()

# inset
pdf("PCA_plots/supp_figure_ComBat_F.pdf",	w=4)
palette(colPal(n))
boxplot(signed*pca3$x[,2]~pd$dates, horiz = T,
	ylab = "",	xaxt="n",ylim = c(-35,65),
	col=1:n)
	
lets = c(letters,LETTERS)[1:n]
text(seq(1,n,by=2), y = 65, lets[seq(1,n,by=2)])
text(seq(2,n,by=2), y = 62, lets[seq(2,n,by=2)])
text(1:n, -33, levels(pd$dates), cex=0.6, srt="270", font=2)
dev.off()


#### 
fitComBat= lmFit(cleanp, mod)
ebComBat = ebayes(fitComBat)
save(fitComBat, ebComBat,file="combat_obj.rda")

## from the other script
fit0 = lmFit(p, mod)
eb0 = ebayes(fit0) # empirical bayes for t-stats
svaobj = sva(p, mod, n.sv=27) # auto is 27
fit1 = lmFit(p, cbind(mod, svaobj$sv))
eb1 = ebayes(fit1)
save(eb1,fit1, file="sva_obj.rda")
###
plot(ebComBat$t[,2], eb1$t[,2])
mean(ebComBat$p[,2] < eb1$p[,2])
plot(-log10(ebComBat$p[,2]), -log10(eb0$p[,2]))
abline(0,1,col="red")
plot(-log10(ebComBat$p[,2]), -log10(eb1$p[,2]))
abline(0,1,col="red")
plot(-log10(eb0$p[,2]), -log10(eb1$p[,2]))
abline(0,1,col="red")

mean(ebComBat$p[,2] < eb1$p[,2])
mean(ebComBat$p[,3] < eb1$p[,3])