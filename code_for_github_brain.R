#### R code for brain preprocessing and normalizing
## andrew jaffe
#

## load packages
library(limma)
library(sva)
library(splines)
library(GEOquery)

# custon functions
# string split wrapper
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), function(y) y[slot])
splitit = function(x) split(seq(along=x),x) # wrapper for split

## plotting function
mypar = function(a=1,b=1,brewer.n=8,brewer.name="Dark2",...){
 par(mar=c(2.5,2.5,1.6,1.1),mgp=c(1.5,.5,0))
 par(mfrow=c(a,b),...)
 palette(brewer.pal(brewer.n,brewer.name))
}

# regress SVs out of expression data
cleaningP = function(y, mod, svaobj,P=ncol(mod)) {
        X=cbind(mod,svaobj$sv)
        Hat=solve(t(X)%*%X)%*%t(X)
        beta=(Hat%*%t(y))
        cleany=y-t(as.matrix(X[,-c(1:P)])%*%beta[-c(1:P),])
        return(cleany)
}

## plot patterns over age
agePlotter = function(y, age, mod, mainText, smoothIt=TRUE,
        orderByAge=TRUE,ylim=NULL,ageBreaks = c(-1, 0, 1, 10, 100),
        ylab="Adjusted Expression",pointColor = 2, lineColor= 1, ...) {

        if(orderByAge) {
                oo = order(age, decreasing = FALSE)
                y = y[oo] ; age = age[oo] ; mod = mod[oo,]
        }

        fit = fitted(lm(y~mod-1))

        fetal = cut(age, breaks =ageBreaks ,lab=F)
        fIndex = splitit(fetal)


        layout(matrix(c(1,1,1,2,2,3,3,4,4,4,4,4),nr = 1,byrow = T))
        palette(brewer.pal(8,"Set1"))

        par(mar = c(4,5,3,0.45))
        if(is.null(ylim)) ylims = range(y,na.rm=TRUE) else ylims = ylim

        xx = jitter(age,amount=0.005)
        plot(y ~ xx,
                subset=fIndex[[1]],
                main = "",      ylab=ylab,xlab="",
                ylim = ylims,cex.axis = 1.5, cex.lab=1.75,
                pch = 21, cex = 1.4,xaxt="n",bg = pointColor,
                xlim=c(range(age[fIndex[[1]]])+c(-0.01,0.07)),...)

        if(smoothIt) lines(age[fIndex[[1]]],fit[fIndex[[1]]],col=lineColor,lwd=6)

        axis(1,at=unique(age[fIndex[[1]]]),
                labels = 40+52*signif(unique(age[fIndex[[1]]]),1), cex.axis=1.5)
        text(x = quantile(age[fIndex[[1]]],0.33), y= min(ylims), "Weeks", cex=1.5)

        # infant + child
        par(mar = c(4, 0.25,3,0.25))
        for(j in 2:3) {
                plot(y ~ age,   subset=fIndex[[j]],
                        main = "",ylab="",xlab="",yaxt = "n", cex=1.4,
                        xlim = range(age[fIndex[[j]]])+c(-0.03,0.03),
                        ylim = ylims, cex.axis = 1.5,pch = 21,  bg=pointColor)

                if(ageBreaks[2] == 0 & smoothIt) lines(age[fIndex[[j]]],fit[fIndex[[j]]],col=lineColor,lwd=6)
                if(ageBreaks[2] < 0 & smoothIt) lines(age[fIndex[[j]]][-1],fit[fIndex[[j]]][-1],col=lineColor,lwd=6)
        }

        # adults
        par(mar = c(4, 0.25,3,1))
        plot(y ~ age,   subset=fIndex[[4]],
                        main = "",ylab="",xlab="",yaxt = "n", cex=1.4,
                        xlim = range(age[fIndex[[4]]])+c(-0.01,0.01),
                        ylim = ylims, cex.axis = 1.5,pch = 21, bg=pointColor)

        if(smoothIt) lines(age[fIndex[[4]]],fit[fIndex[[4]]],col=lineColor,lwd=6)

        mtext(mainText, outer=T, line=-2.5,cex=1.35)

        mtext("Age", side=1, outer=T, line=-1.5,cex=1.35)
}

###################################
### PREPROCESS FROM GEO
# to get phenotype and processed data
theData = getGEO("GSE30272")

pdList = lapply(theData,pData)
pd = do.call("rbind",pdList)
pd = pd[,c(1,2,grep("characteristics", names(pd)))]
pd = pd[,1:10]

pList = lapply(theData, exprs)
p = do.call("cbind", pList)

colnames(p) = rownames(pd) = pd$title

# do some formatting
colNames = apply(pd[,3:10], 2, function(x) {
	ss(as.character(x),": ", 1)[1]
})
colnames(pd)[3:10] = colNames

pd[,3:10] = apply(pd[,3:10], 2, function(x) {
	ss(as.character(x),": ", 2)
})
pd$title = as.character(pd$title)
pd$geo_accession = as.character(pd$geo_accession)
for(i in c(3,4,7,8,9)) pd[,i] = as.numeric(pd[,i])

## map
map = as(featureData(theData[[1]]), "data.frame")
rownames(map) = rownames(p) = map$OligoID

## RAW
getGEOSuppFiles("GSE30272")
system("gunzip GSE30272/GSE30272_RGna.n269.o30176.log2ratio.loess.MADout.KNNimp.txt.gz")
p = read.table("GSE30272/GSE30272_RGna.n269.o30176.log2ratio.loess.MADout.KNNimp.txt",
	header=TRUE, as.is=TRUE, row.names=1)
	
map = map[rownames(p),]

## order by age
oo = order(pd$age)
p = p[,oo] ; pd = pd[oo,]
p = as.matrix(p)

######################
##### perform modeling
# fetal offset
fetal = ifelse(pd$age < 0, 1, 0)

# linear spline
age = pd$age
ageF = age
ageF[ageF <0]=0

modLin = model.matrix(~age+ageF+fetal)
svaobjLin = sva(p, modLin)

## spline, deg=2
ageSpline = bs(pd$age, knots = 0, degree=2)
modCurve = model.matrix(~ageSpline + fetal)
svaobjCurve = sva(p, modCurve)

ageSpline1 = bs(pd$age[pd$age < 0], degree=2)
ageSpline1 = cbind(ageSpline1, 0,0,0,0,0,0)
ageSpline2 = cbind(0,0,bs(pd$age[pd$age > 0], knots = c(0.6,10,20,50),degree=2))
mod = model.matrix(~rbind(ageSpline1, ageSpline2)+fetal)
svaobj = sva(p, mod)

save(svaobjLin,svaobjCurve,svaobj, 
	file="svaobj_fetalBreaks_svaPaper.rda")

# load("svaobj_fetalBreaks_svaPaper.rda")

#########################################
# check residuals from taking out each model from raw data

resLin = p - fitted(lmFit(p,modLin))
resCurve = p - fitted(lmFit(p,modCurve))
resRaw = p - rowMeans(p)
resMod = p - fitted(lmFit(p,mod))

rssLin = rowSums(resLin^2)
rssCurve = rowSums(resCurve^2)
rssRaw = rowSums(resRaw^2)
rssMod = rowSums(resMod^2)

infantChild = which(pd$age > 0 & pd$age < 20)
rssLinIC = rowSums(resLin[,infantChild]^2)
rssCurveIC = rowSums(resCurve[,infantChild]^2)
rssRawIC = rowSums(resRaw[,infantChild]^2)
rssModIC = rowSums(resMod[,infantChild]^2)

######
cleanpLin = cleaningP(p,modLin,svaobjLin)
cleanpCurve = cleaningP(p,modCurve,svaobjCurve)
cleanp = cleaningP(p,mod,svaobj)


### most "spliny"
oo = order(rssMod - rssCurve, decreasing=FALSE)
# age = log2(pd$age+1)
i = oo[1]

pdf("plots/brainSpline_rawLine.pdf", h = 4.5, w = 6)
agePlotter(p[i,], pd$age, modLin, mainText=rownames(p)[i],ylim=range(p[i,]))
dev.off()

pdf("plots/brainSpline_rawCurve.pdf", h = 4.5, w = 6)
agePlotter(p[i,], pd$age, modCurve, mainText=rownames(p)[i],ylim=range(p[i,]))
dev.off()

pdf("plots/brainSpline_rawFinal.pdf", h =  4.5, w = 6)
agePlotter(p[i,], pd$age, mod, mainText=rownames(p)[i],ylim=range(p[i,]))
dev.off()

pdf("plots/brainSpline_cleanLine.pdf", h = 4.5, w = 6)
agePlotter(cleanpLin[i,], pd$age, mod, mainText=rownames(p)[i],
	 ylim=range(p[i,]), lineColor = "green")
dev.off()

pdf("plots/brainSpline_cleanCurve.pdf", h = 4.5, w = 6)
agePlotter(cleanpCurve[i,], pd$age, mod, mainText=rownames(p)[i],
	 ylim=range(p[i,]), lineColor = "green")
dev.off()

pdf("plots/brainSpline_cleanFinal.pdf", h = 4.5, w = 6)
agePlotter(cleanp[i,], pd$age, mod, mainText=rownames(p)[i],
	ylim=range(p[i,]), lineColor = "green")
dev.off()


########
## look at residuals
cleanpDiff = cleanp - p
cleanpLinDiff = cleanpLin - p
cleanpCurveDiff = cleanpCurve - p

rss2 = rowSums(cleanpDiff^2)
rssLin2 = rowSums(cleanpLinDiff^2)
rssCurve2 = rowSums(cleanpCurveDiff^2)

mypar(brewer.name="Set1")
plot(density(log10(rss2)), lwd = 3, col = 1)
lines(density(log10(rssLin2)), lwd = 3, col = 2)
lines(density(log10(rssCurve2)), lwd = 3, col = 3)

smoothScatter(log10(rss2), log10(rssLin2),
	xlab = "Final",ylab="Line");abline(0,1,col="red")
smoothScatter(log10(rss2), log10(rssCurve2),
	xlab = "Final",ylab="Curve");abline(0,1,col="red")
smoothScatter(log10(rssLin2), log10(rssCurve2),
	xlab = "Line",ylab="Curve");abline(0,1,col="red")

pdf("smooth_scatters_brainModels_supplementary.pdf")
M = log10(rss2)  - log10(rssLin2)
A = (log10(rss2)  + log10(rssLin2))/2
smoothScatter(M~A,	ylab = "Log10(Final) - Log10(Line)",
	xlab="Avg Log10(RSS)",	ylim = c(-0.8, 0.45),
	cex.axis=1.6)
abline(h=mean(M), col="red",lwd=1.5)
abline(h=0, col="black",lty=2,lwd=1.5)
text(x = 3, y = 0.1, sum(M > 0),cex=2)
text(x = 3, y =- 0.1, sum(M < 0),cex=2)


M = (log10(rss2)  - log10(rssCurve2))
A = (log10(rss2)  + log10(rssLin2))/2
smoothScatter(M~A,	ylab = "Log10(Final) - Log10(Curve)",
	xlab="Avg Log10(RSS)",	ylim = c(-0.8, 0.45),
	cex.axis=1.6)
abline(h=mean(M), col="red",lwd=1.5)
abline(h=0, col="black",lty=2,lwd=1.5)
text(x = 3, y = 0.1, sum(M > 0),cex=2)
text(x = 3, y =- 0.1, sum(M < 0),cex=2)

M = (log10(rssCurve2)  - log10(rssLin2)) 
A = (log10(rss2)  + log10(rssLin2))/2
smoothScatter(M ~ A, ylab = "Log10(Curve) - Log10(Line)",
	xlab="Avg Log10(RSS)",	ylim = c(-0.8, 0.45),
	cex.axis=1.6)
abline(h=mean(M), col="red",lwd=1.5)
abline(h=0, col="black",lty=2,lwd=1.5)
text(x = 3, y = 0.1, sum(M > 0),cex=2)
text(x = 3, y =- 0.1, sum(M < 0),cex=2)
dev.off()

###### PCA
pcaRaw = prcomp(t(p))
pca = prcomp(t(cleanp))
pcaLin = prcomp(t(cleanpLin))
pcaCurve = prcomp(t(cleanpCurve))

pcaVarsRaw=signif(((pcaRaw$sdev)^2)/(sum((pcaRaw$sdev)^2)),3)*100
pcaVars=signif(((pca$sdev)^2)/(sum((pca$sdev)^2)),3)*100
pcaVarsLin=signif(((pcaLin$sdev)^2)/(sum((pcaLin$sdev)^2)),3)*100
pcaVarsCurve=signif(((pcaCurve$sdev)^2)/(sum((pcaCurve$sdev)^2)),3)*100



ageGroup2 = cut(pd$age, breaks = c(-0.5,0,0.6,seq(10,80,by=10)), labels=F)
ageGroup2Names = c("Fetal","Infant","Child","Teens","20s","30s",
	"40s","50s","60s","70+s")
palette(c("purple","darkblue","blue","lightblue","grey","green",
	"yellow","orange","darkorange","red"))

pdf("pca_by_diff_models.pdf", h = 12, w = 8)
mypar(3,2)
palette(c("purple","darkblue","blue","lightblue","grey","green",
	"yellow","orange","darkorange","red"))

plot(pca$x[,1], pca$x[,2], pch = 21, bg = as.numeric(ageGroup2),
	main = "Final Spline Model",cex=1.2,
	xlab = paste("PC1:",pcaVars[1],"% of Variance Explained"),
	ylab = paste("PC2:",pcaVars[2],"% of Variance Explained"))
plot(pca$x[,3], pca$x[,2], pch = 21, bg = as.numeric(ageGroup2),
	main = "Final Spline Model",cex=1.2,
	ylab = paste("PC2:",pcaVars[2],"% of Variance Explained"),
	xlab = paste("PC3:",pcaVars[3],"% of Variance Explained"))

plot(pcaLin$x[,1], pcaLin$x[,2], pch = 21, bg = as.numeric(ageGroup2),
	main = "Linear Spline Model",cex=1.2,
	xlab = paste("PC1:",pcaVarsLin[1],"% of Variance Explained"),
	ylab = paste("PC2:",pcaVarsLin[2],"% of Variance Explained"))
plot(pcaLin$x[,3], pcaLin$x[,2], pch = 21, bg = as.numeric(ageGroup2),
	main = "Linear Spline Model",cex=1.2,
	ylab = paste("PC2:",pcaVarsLin[2],"% of Variance Explained"),
	xlab = paste("PC3:",pcaVarsLin[3],"% of Variance Explained"))

plot(pcaCurve$x[,1], pcaCurve$x[,2], pch = 21, bg = as.numeric(ageGroup2),
	main = "Curve Spline Model",cex=1.2,
	xlab = paste("PC1:",pcaVarsCurve[1],"% of Variance Explained"),
	ylab = paste("PC2:",pcaVarsCurve[2],"% of Variance Explained"))
plot(pcaCurve$x[,3], pcaCurve$x[,2], pch = 21, bg = as.numeric(ageGroup2),
	main = "Curve Spline Model",cex=1.2,
	ylab = paste("PC2:",pcaVarsCurve[2],"% of Variance Explained"),
	xlab = paste("PC3:",pcaVarsCurve[3],"% of Variance Explained"))
dev.off()

pdf("pcs_vs_age_by_diff_models_Final.pdf", h = 5, w = 8)
for(i in 1:6) agePlotter(pca$x[,i], pd$age, mod,
	mainText=paste("Final Model, PC",i), ylab=paste0("PC",i))
dev.off()

pdf("pcs_vs_age_by_diff_models_Raw.pdf", h = 5, w = 8)
for(i in 1:6) agePlotter(pcaRaw$x[,i], pd$age, mod, 
	mainText=paste("Raw, PC",i), ylab=paste0("PC",i))
dev.off()

pdf("pcs_vs_age_by_diff_models_Line.pdf", h = 5, w = 8)
for(i in 1:6) agePlotter(pcaLin$x[,i], pd$age, mod,
	mainText=paste("Linear Spline Model, PC",i), ylab=paste0("PC",i))
dev.off()

pdf("pcs_vs_age_by_diff_models_Curve.pdf", h = 5, w = 8)
for(i in 1:6) agePlotter(pcaCurve$x[,i], pd$age, mod, 
	mainText=paste("Curve Spline Model, PC",i), ylab=paste0("PC",i))
dev.off()