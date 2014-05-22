#########
# Preprocess StemCellDB in GEO
# Andrew Jaffe
# Updated 4/8/12
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
