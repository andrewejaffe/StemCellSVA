#

library(GEOquery)
ss = function(x, pattern, slot=1,...) sapply(strsplit(x,pattern,...), function(y) y[slot])

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

###### THIS IS WHAT SHOULD WORK, ABOVE


## RAW
getGEOSuppFiles("GSE30272")
# system("gunzip GSE30272/GSE30272_RGna.n269.o30176.log2ratio.loess.MADout.KNNimp.txt.gz")
p = read.table("GSE30272/GSE30272_RGna.n269.o30176.log2ratio.loess.MADout.KNNimp.txt",
	header=TRUE, as.is=TRUE, row.names=1)
	
map = map[rownames(p),]

## order by age
oo = order(pd$age)
p = p[,oo] ; pd = pd[oo,]

save(p,map,pd, file="brain_data_from_geo.rda")

