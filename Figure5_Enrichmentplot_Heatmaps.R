############################################################################################################################################################
#Figure 5: Reanalysis of Genome biology data for genes ranked by gene expression differences between current and former smokers
#Entrez V11.1
#
############################################################################################################################################################
ncfile <- read.delim("Paper/Files/GSEA/GSE7895/Gen_biol_rmaexpfile1_gs.txt")
ncDemo <- read.delim("Paper/Files/GSEA/GSE7895/Genbiol.txt")
nc.genes<- rownames(ncfile)

ncpatient <- ncDemo[,1]

pat<-c()
for(i in 1:length(ncpatient)){													###arrange exp and Demo file in same order
indi<-grep(ncpatient[i],colnames(ncfile))
pat<- append(pat,indi)}
#Arrange Exp accorinf to indices
ncfile2 <- ncfile[,pat]
ncDemo1<-ncDemo


f.ind<-which(as.character(ncDemo1[,7])=="Former Smoker") #31
quitt<-ncDemo1[f.ind,5]
quitt1<-as.numeric(as.character( quitt ))

fExp<-ncfile2[,f.ind]
f <- rep(1,length(f.ind))
formExp<-fExp


#current smokencrs

c.ind<-which(as.character(ncDemo1[,7])=="Current Smoker") #52
c <- rep(0,length(c.ind))
currExp<-ncfile2[,c.ind]


mat <- cbind(currExp,formExp)



ncDemo2 <- ncDemo1[c(c.ind,f.ind),]
allPatients <- c(colnames(currExp),colnames(formExp))
age  <-ncDemo2[,4]
status <- c(c,f)
res <-c()
library(nlme)

for(i in 1:nrow(mat)){

gene.exp<-as.numeric(mat[i,])
try({

m1<-lm(gene.exp~as.factor(status) + as.numeric(age))


#calculate and see the p-value for factor y
coeff<-summary(m1)$coefficients[2,3]

pval<-anova(m1)[1,5]		

	

cp<-cbind(as.character(nc.genes[i]),coeff,pval)
print (i)
res<-rbind(res,cp)
},TRUE)
rownames(res)<-c()
}


rowname <- res[,1]
res1 <- res[,2:3]
rownames(res1) <- rowname
#write down file
write.table(res1,"Paper/Files/GSEA/GSE7895/GSE7895_CvsF_model_results.txt",sep="\t")
write.table(res1[,1:2],"Paper/Files/GSEA/GSE7895/072912_GB_CvsF_model_ranklist.rnk",sep="\t")

############################################################################################################################################################
#Figure 5A: Enrichment plot 
#
############################################################################################################################################################
source("gsea.R")
pdf ("Figure5A: GSE7895_CvsFRank_SC_geneset_enrichmentplot.pdf")
plot.gsea("Paper/Files/GSEA/GSE7895/072912_GB_CvsF_model_ranklist.rnk","Paper/Files/GSEA/GSE7895/081312_GSEA_GB_CvsFrank_SC_geneset_enrichmentplotfile.txt")
dev.off()

############################################################################################################################################################
#Figure 5B,C: Heatmaps
#
############################################################################################################################################################
