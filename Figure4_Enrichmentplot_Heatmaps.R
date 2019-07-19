############################################################################################################################################################
#Figure 4: Reanalysis of Genome biology data for genes ranked by gene expression differences between current and former smokers
#Entrez V11.1
#


############################################################################################################################################################
ncfile <- read.delim("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/GSEA/GSE7895/Gen_biol_rmaexpfile1_gs.txt")
ncDemo <- read.delim("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/GSEA/GSE7895/Genbiol.txt")
samples<- colnames(ncfile)
samples <- samples[2:105]
nc.genes<- ncfile[,2]
ncfile <-ncfile[,3:106]
ncpatient <- ncDemo[,1]
colnames(ncfile)<-c()
colnames(ncfile) <- samples
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
f <- rep(0,length(f.ind))
formExp<-fExp


#current smokencrs

c.ind<-which(as.character(ncDemo1[,7])=="Current Smoker") #52
c <- rep(1,length(c.ind))
currExp<-ncfile2[,c.ind]

#Never smokers
n.ind<-which(as.character(ncDemo1[,7])=="Never Smoker") 
n <- rep(0,length(n.ind))
nevExp<-ncfile2[,n.ind]


mat <- cbind(currExp,formExp)
mat1<- cbind(currExp,formExp,nevExp)



##############################################Reanalysis of Genome biology for GSEA ranklist#########################################
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
write.table(res1,"Paper/Files/GSEA/GSE7895/GSE7895_CvsF_10_model_results.txt",sep="\t")
write.table(res1[,1:2],"Paper/Files/GSEA/GSE7895/072912_GB_CvsF_model_ranklist.rnk",sep="\t")
#####################Figure 4B Heatmap###################################################################################

setwd("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/")
lead.file <- read.table("GSEA/GSE7895/GSEA_SC_119genes_genomeBio_bronch_ranklist_LeadingEdgeGenes.txt",sep="\t")
genes.lead<-lead.file[,1] #30 genes
ind.gb.genes<- match(as.character(genes.lead),nc.genes,nomatch=0)

mat2 <- mat1[ind.gb.genes,]
nc.genes[ind.gb.genes]
rownames(mat2)<- nc.genes[ind.gb.genes]
library(heatmap3)
#Load graphics/plotting packages that has red/blue colors for the heat map
library(gplots)
#create a variable with colors that will be used in heat map
bluered<-bluered(256)
##Creating a vector that contains colorsrepresenting the class of each sample
x<-c("darkblue","darkgreen","yellow")
#Tclasslables<-c(rep("red",vsamples[1]),rep("blue",vsamples[2]),rep("green",vsamples[3]),rep("purple",vsamples[4]),rep("orange",vsamples[5]))
Tclasslables<-c(rep(x[1],length(c.ind)),rep(x[2],length(f.ind)),rep(x[3],length(n.ind)))


pdf(file="/protected/projects/pulmarray/U01_SmokingCessation/Paper/Figures/Figure4C LeadingEdge_Heatmap_GenomeBiology_with_SC.pdf")	
heatmap3(mat2, main="", col=bluered,

         row.clustering="supervised", col.clustering="supervised",

         row.dendrogram=TRUE, ColSideColors=Tclasslables)

#close output file
dev.off()	



############################################################################################################################################################
#Figure 4A: Enrichment plot 
#
############################################################################################################################################################
source("/protected/projects/pulmarray/U01_SmokingCessation/Paper/R_Programs_Paper/gsea.R")
pdf ("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Figures/Figure 4_GSE7895_GenomeBio_CvsFRank_1_vs_0_SC_geneset_enrichmentplot.pdf")
plot.gsea("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/GSEA/GSE7895/GSE7895_CvsF_10_model_ranklist.rnk", "/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/GSEA/GSE7895/081013_GSEA_GB_CvsFrank_10_SC_geneset_enrichmentplotfile.txt")
dev.off()

############################################################################################################################################################
#Figure 4B,C: Leading edge Heatmaps
#
############################################################################################################################################################
setwd("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/ModelResults")
Exp<-read.delim("SC_ExpressionData_All_genes_genesymbol.txt")
genes <- rownames(Exp)
dim(Exp)
setwd("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/")
DemoBmc <- read.table ("BMC_SmokingCessation_Demographics.txt",sep="\t",header=TRUE)
DemoMn <- read.table ("MN_SmokingCessation_Demographics.txt",sep="\t",header=TRUE)
dim(DemoBmc)
#25 13
dim(DemoMn)
#[1] 32 13

###########read leading edge genes
lead.file <- read.table("GSEA/GSE7895/GSEA_SC_119genes_genomeBio_bronch_ranklist_LeadingEdgeGenes.txt",sep="\t")
genes.lead<-lead.file[,1] #30 genes


#Create vector of patient ids, excluding outliers (N208_24)
#Create vector of patient ids
patient.bmc <- c(107,116,206)
patient.mn <- c(205,208,227,243,265)


#####Arrange Demo files BMC##########################for excluding patients
patient.ind.bmc=c()
for(i in 1:length(patient.bmc)){
bmc<-grep(patient.bmc[i],DemoBmc[,1])
patient.ind.bmc <- append(patient.ind.bmc,bmc)}
DemoBmc1 <- DemoBmc[patient.ind.bmc,]

####Arrange Demo file Mn
patient.ind.bmc=c()
for(i in 1:length(patient.mn)){
mn<-grep(patient.mn[i],DemoMn[,2])
patient.ind.bmc <- append(patient.ind.bmc,mn)}
DemoMn1 <- DemoMn[patient.ind.bmc,]

cell.id<-DemoMn1[,1]
patient <- c(patient.bmc,cell.id)


#Extract out colnames that match patient vector for ordering of matrice
patient.ind.bmc=c()
for(i in 1:length(patient)){
ind.bmc<-grep(patient[i],colnames(Exp))
patient.ind.bmc <- append(patient.ind.bmc,ind.bmc)}
#Arrange Exp accorinf to indices of demo files
Exp2 <- Exp[,patient.ind.bmc]


###########################################################end of excd patients
####Create vector corresponding to number of time points per patient

All_samples <- c(as.character(DemoBmc1[,1]),as.character(DemoMn1[,2]))

#excl_samples <- c("N205_4"," N205_8  ")
excl_samples <- c("N205_4"," N205_8  ","N227_16")

ex.ind <- match(excl_samples,All_samples,nomatch=0)
Mn.ind <-match(excl_samples,as.character(DemoMn1[,2]),nomatch=0)
bmc.ind<- match(excl_samples,as.character(DemoBmc1[,1]),nomatch=0)
###Exclude samples

B<-DemoBmc1
D <- DemoMn1[-Mn.ind,]
DemoMn2 <-D

DemoBmc2<-B

vsamples<- c(2,3,5,3,5,4,5,4)
####Make vector of patient id's
patient_id <- c(DemoBmc2[,5],DemoMn2[,3])
patient_id <-as.factor(patient_id)
Time <- c(DemoBmc2[,3],DemoMn2[,4]) # 
RIN<-c(DemoBmc2[,13],DemoMn2[,9])
RIN<-as.numeric(RIN)


colnames(Exp)<-c()
names<-c(as.character(DemoBmc2[,1]),as.character(DemoMn2[,2]))
colnames(Exp) <- names



#############################################reduce nasal matrice to leading edge genes

ind.genes<- match(as.character(genes.lead),genes,nomatch=0)

NewExp <- Exp[c(ind.genes),]

Matrix<-as.matrix(NewExp)														  #Adjust residuals for patient effect and RIN

name<-colnames(Matrix)
model <- model.matrix(~1 + patient_id)										  #Get residuals after adjusting for patient and RIN

fileGE<-as.matrix(Matrix) %*% (diag(dim(Matrix)[2]) - model %*% solve(t(model) %*% model) %*% t(model)) #Get residuals after adjusting for patient and RIN

colnames(fileGE) <- colnames(Matrix)

Time <- c(DemoBmc2[,3],DemoMn2[,4])
time.ind <-order(Time)

fileGE1<-fileGE[,time.ind]


library(heatmap3)																			#Load graphics/plotting packages that has red/blue colors for the heat map

library(gplots)																				#Create a variable with colors that will be used in heat map

bluered<-bluered(256)                                                                       ##Creating a vector that contains colorsrepresenting the class of each sample
index <- length(fileGE[,1])

vsamples<-c(8,6,7,6,6)                                                                      #Vector representing number of samples per time
P1<-c()																						#Add white lines between time points on heatmap
P2<-c()
P3<-c()
P4<-c()
P5<-c()
P1<- fileGE1[,1:8]
P1<-cbind(P1,c(rep(0,length(index))))
P2<- fileGE1[,9:14]
P2<-cbind(P2,c(rep(0,length(index))))
P3<- fileGE1[,15:21]
P3<-cbind(P3,c(rep(0,length(index))))
P4<- fileGE1[,22:27]
P4<-cbind(P4,c(rep(0,length(index))))
P5<- fileGE1[,28:33]

E<-cbind(P1,P2,P3,P4,P5)
#################################################Heatmap###################################################
vsamples<-c(8,1,6,1,7,1,6,1,6)
#vsamples<-c(8,6,7,6,6)
library(heatmap3)
#Load graphics/plotting packages that has red/blue colors for the heat map
library(gplots)
#create a variable with colors that will be used in heat map
bluered<-bluered(256)
##Creating a vector that contains colorsrepresenting the class of each sample
		Tclasslables<-c(rep("darkblue",vsamples[1]),rep("white",vsamples[2]),rep("blue",vsamples[3]),rep("white",vsamples[4]),
		rep("darkgreen",vsamples[5]),rep("white",vsamples[6]),rep("pale green",vsamples[7]),rep("white",vsamples[8]),rep("yellow",vsamples[9]))

pdf(file="/protected/projects/pulmarray/U01_SmokingCessation/Paper/Figures/Figure4B LeadingEdge_Heatmap_SmokingCessation_with_GB.pdf")	
heatmap3(E, main="", col=bluered,

         row.clustering="supervised", col.clustering="supervised",

         row.dendrogram=TRUE, ColSideColors=Tclasslables)

#close output file
dev.off()	




