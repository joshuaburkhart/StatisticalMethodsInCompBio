#A. Change working directory to where the lab files are located. Depending on what you have installed, you may need to add additional packages here.

if(FALSE){
source("http://www.bioconductor.org/biocLite.R")
biocLite("siggenes")
biocLite("RankProd")
biocLite("limma")
biocLite("fibroEset")
biocLite("made4")
biocLite("annaffy")
biocLite("hgu95av2.db")
}

require(affy)
require(annaffy)
require(hgu95av2.db)
require(made4)


# B. Load data and examine annotation 
data.vsn<- read.csv("/Users/joshuaburkhart/SoftwareProjects/StatisticalMethodsInCompBio/InClass/data.vsn.csv", as.is=TRUE, row.names=1)
dim(data.vsn)
annt<-read.table("/Users/joshuaburkhart/SoftwareProjects/StatisticalMethodsInCompBio/InClass/annt.txt", header=TRUE)
annt[1:2,]

annt$Donor
table(annt$Donor)
table(annt$Gender)

names(data.vsn)
annt
rownames(annt) <-annt$Cels

#C. Create Bioconductor eSet object from files you loaded 

makeEset<-function(eSet, annt){
    #Creating an ExpressionSet from eSet, a normalized gene expression matrix
    # and annt, a data.frame containing annotation
    metadata <- data.frame(labelDescription = colnames(annt), row.names=colnames(annt))
    phenoData<-new("AnnotatedDataFrame", data=annt, varMetadata=metadata)
    if (inherits(eSet, "data.frame")) eSet= as.matrix(eSet)
    if (inherits(eSet, "ExpressionSet")) eSet=exprs(eSet)
    data.eSet<-new("ExpressionSet", exprs=eSet, phenoData=phenoData)
    print(varLabels(data.eSet))
    return(data.eSet)
}

eSet<-makeEset(data.vsn, annt)

#D. We will set this up for a human (HSA) vs non-human comparison 
human<- eSet$Donor=="Hsa"
table(human)
eSet$Human<-human
eSet


#E. Retrieving Gene Annotation for array 
affy.id = featureNames(eSet)
affy.symbols<-aafSymbol(affy.id, "hgu95av2.db")
affy.symbols <-getText(affy.symbols)
names(affy.symbols)<-featureNames(eSet)


# F. Using Limma for linear models per gene/transcript

require(limma)

design= model.matrix(~eSet$Human)
fit <- lmFit(eSet,design)
fit <- eBayes(fit)
topTable(fit,coef=2)

limmaRes = topTable(fit,coef=2,  p.value=0.001, number=500)
print(nrow(limmaRes))

# G. Visualizing the Limma results

heatplot(eSet[rownames(limmaRes),], classvec=eSet$Donor, labRow=affy.symbols[limmaRes$ID], labCol=eSet$Donor)

par(mfrow=c(2,1))
qqt(fit$t[,2],df=fit$df.residual+fit$df.prior)
abline(0,1)
volcanoplot(fit,coef=2,highlight=2)


# H. RankProduct 
require(RankProd)
RP.out <- RP(eSet, eSet$Human, rand=123)
plotRP(RP.out, cutoff=0.05)
RP.res = topGene(RP.out,cutoff=0.05,method="pfp",logged=TRUE,logbase=2,gene.names=affy.symbols)
names(RP.res)
RP.res$Table1[1:10,]
RP.res$Table2[1:10,]


#I: Compare to full data set
require(fibroEset)
data(fibroEset)
phenoData(fibroEset)






