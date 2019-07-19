#######################Figure 3A Enrichment plot GSE16008################################
source("/protected/projects/pulmarray/U01_SmokingCessation/Paper/R_Programs_Paper/gsea.R")
pdf ("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Figures/Figure 3A_GSE16008_PhysiolGen_CvsNMainRank_SC_geneset_enrichmentplot.pdf")
plot.gsea("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/GSEA/GSE16008/PhysiGenomics..13vs14.Nose.17881g.2batch.random.sort.rnk","/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/GSEA/GSE16008/NoseEnrichment.txt")
dev.off()



##Figure 3B,C Leading Edge Boxplots across two datasets: SmokingCessation and GSE16008

setwd("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/")
Exp<-read.delim("ModelResults/SC_ExpressionData_All_genes_genesymbol.txt")
genes <- rownames(Exp)
dim(Exp) #19722 33

###Read Demographics

DemoBmc <- read.table ("BMC_SmokingCessation_Demographics.txt",sep="\t",header=TRUE)
DemoMn <- read.table ("MN_SmokingCessation_Demographics.txt",sep="\t",header=TRUE)
dim(DemoBmc)
#25 13
dim(DemoMn)
#[1] 32 13
lead.file <- read.table("GSEA/GSE16008/GSEA_SC_119genes_PhysiolGenom_nose_ranklist_LeadingEdgeGenes.txt",sep="\t")
genes.lead<-lead.file[,1] #10 genes


patient.bmc <- c(107,116,206)																	###Create vectors of patient ids to include in analysis
patient.mn <- c(205,208,227,243,265)


                                                                                                ###Arrange Demo files BMC
patient.ind.bmc=c()
for(i in 1:length(patient.bmc)){
bmc<-grep(patient.bmc[i],DemoBmc[,1])
patient.ind.bmc <- append(patient.ind.bmc,bmc)}
DemoBmc1 <- DemoBmc[patient.ind.bmc,]

																							   ###Arrange Demo file Mn
patient.ind.bmc=c()
for(i in 1:length(patient.mn)){
mn<-grep(patient.mn[i],DemoMn[,2])
patient.ind.bmc <- append(patient.ind.bmc,mn)}
DemoMn1 <- DemoMn[patient.ind.bmc,]

cell.id<-DemoMn1[,1]
patient <- c(patient.bmc,cell.id)


																					     ###Extract out colnames that match patient vector for ordering of matrice
																						 ###Create vector corresponding to number of time points per patient

All_samples <- c(as.character(DemoBmc1[,1]),as.character(DemoMn1[,2]))
excl_samples <- c("N205_4"," N205_8  ","N227_16")											   ###Exclude outliars
ex.ind <- match(excl_samples,All_samples,nomatch=0)
Mn.ind <-match(excl_samples,as.character(DemoMn1[,2]),nomatch=0)
bmc.ind<- match(excl_samples,as.character(DemoBmc1[,1]),nomatch=0)
B<-DemoBmc1
D <- DemoMn1[-Mn.ind,]
DemoMn2 <-D
DemoBmc2<-B
vsamples<- c(2,3,5,3,5,4,5,4)
																								###Make vector of patient id's
patient_id <- c(DemoBmc2[,5],DemoMn2[,3])
patient_id <-as.factor(patient_id)
Time <- c(DemoBmc2[,3],DemoMn2[,4]) # 
RIN<-c(DemoBmc2[,13],DemoMn2[,9])
RIN<-as.numeric(RIN)
colnames(Exp)<-c()
names<-c(as.character(DemoBmc2[,1]),as.character(DemoMn2[,2]))
colnames(Exp) <- names




match.ind<- match(as.character(genes.lead),genes,nomatch=0)									### match genes

fileGE <- Exp[match.ind,]
rownames(fileGE)<- genes.lead																            
														                ### vector for time points


Matrix<-as.matrix(fileGE)														  #Adjust residuals for patient effect and RIN

name<-colnames(Matrix)
model <- model.matrix(~1 + patient_id)										  #Get residuals after adjusting for patient and RIN

fileGE1<-as.matrix(Matrix) %*% (diag(dim(Matrix)[2]) - model %*% solve(t(model) %*% model) %*% t(model)) #Get residuals after adjusting for patient and RIN

colnames(fileGE1) <- colnames(Matrix)

Time <- c(DemoBmc2[,3],DemoMn2[,4])
time.ind <-order(Time)

fileGE2<-fileGE1[,time.ind]

Time1<- Time[time.ind]

vsamples<-c(8,6,7,6,6)	

library(heatmap3)																			#Load graphics/plotting packages that has red/blue colors for the heat map

library(gplots)																				#Create a variable with colors that will be used in heat map

bluered<-bluered(256)                                                                       ##Creating a vector that contains colorsrepresenting the class of each sample
index <- length(fileGE2[,1])

vsamples<-c(8,6,7,6,6)                                                                      #Vector representing number of samples per time
P1<-c()																						#Add white lines between time points on heatmap
P2<-c()
P3<-c()
P4<-c()
P5<-c()
P1<- fileGE2[,1:8]
P1<-cbind(P1,c(rep(0,length(index))))
P2<- fileGE2[,9:14]
P2<-cbind(P2,c(rep(0,length(index))))
P3<- fileGE2[,15:21]
P3<-cbind(P3,c(rep(0,length(index))))
P4<- fileGE2[,22:27]
P4<-cbind(P4,c(rep(0,length(index))))
P5<- fileGE2[,28:33]

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

pdf(file="/protected/projects/pulmarray/U01_SmokingCessation/Paper/Figures/Figure3B LeadingEdge_Nose_Heatmap_SmokingCessation_with_GSE16008.pdf")	
heatmap3(E, main="", col=bluered,

         row.clustering="unsupervised", col.clustering="supervised",

         row.dendrogram=TRUE, ColSideColors=Tclasslables)

#close output file
dev.off()	





##################GSE16008##############


#read in Physiol Genomics data
Exp.p <- read.delim("GSEA/GSE16008/Aim01_physiol_demographics_CvsN_combined_expressionData.txt")
trans.exp <- Exp.p[,1]
Exp <- Exp.p[,2:53]
Demo.p <- read.delim("GSEA/GSE16008/Aim01_physiol_demographics_CvsN.txt")
#Demo <- read.delim("Aim01_physiol_demographics_CvsN_bronch.txt")
conv.data <- read.table("/protected/individuals/kahkshan/U01-aim1/nose/TranscriptIds_to_genesymbols_aim01_dashremoved.txt")
gene.conv <- unique(conv.data[,2])
u.ind <- match(gene.conv,conv.data[,2],nomatch=0)
conv.data1 <- conv.data[u.ind,]
c.ind <- match(trans.exp,conv.data1[,1],nomatch=0)
conv.data2 <-conv.data1[c.ind,]
genes.p <- conv.data2[,2]
exp.ind <- match(conv.data2[,1],trans.exp)
Exp1.p<- Exp[exp.ind,]
rownames(Exp1.p) <- genes.p
batch <- as.factor(Demo.p$Batch)
###residuals after adjusting for batch
Matrix<-as.matrix(Exp1.p)														  #Adjust residuals for patient effect and RIN

name<-colnames(Matrix)
model <- model.matrix(~1 + batch)									  #Get residuals after adjusting for patient and RIN

fileGE<-as.matrix(Matrix) %*% (diag(dim(Matrix)[2]) - model %*% solve(t(model) %*% model) %*% t(model)) #Get residuals after adjusting for patient and RIN

colnames(fileGE) <- colnames(Matrix)



c.ind.p<-which(as.character(Demo.p$Smokest)=="Current") #
n.ind.p<-which(as.character(Demo.p$Smokest)=="Never") #
b.ind <- which(as.character(Demo.p$site)=="Bronch")
n.ind <- which(as.character(Demo.p$site)=="Nose")


curr.bronch.ind <- intersect(c.ind.p,b.ind)
never.bronch.ind <- intersect(n.ind.p,b.ind)

curr.nose.ind <- intersect(c.ind.p,n.ind)
never.nose.ind <- intersect(n.ind.p,n.ind)

g.ind.p <- match(genes.lead,rownames(Exp1.p),nomatch=0)
Exp2.p<- Exp1.p[g.ind.p,c(curr.nose.ind,never.nose.ind)]#27:53

## Run prediction code
# source("/protected/projects/pulmarray/U01_SmokingCessation/Paper/R_Programs_Paper/svd_project.R")
# res = svd.project(train.data=fileGE1, test.data=Exp3.p)

## Get scores
# train.scores = res$train.scores
# test.scores = res$test.scores
###Heatmap GSE16008

library(heatmap3)                                            #Load graphics/plotting packages that has red/blue colors for the heat map

library(gplots)                                              #create a variable with colors that will be used in heat map

bluered<-bluered(256)
source("/protected/projects/pulmarray/U01_SmokingCessation/Paper/R_Programs_Paper/zscore_rows.R")
#new<-zscore.rows(Exp2.p)
gluc.class <- c(rep("Darkblue",length(curr.nose.ind)), rep("yellow",length(never.nose.ind)))#######Zscore the new matrice@############

pdf(file="/protected/projects/pulmarray/U01_SmokingCessation/Paper/Figures/Figure3C Heatmap of leading edge genes SmokingCessation Nose GSE16008.pdf")

heatmap3(Exp2.p, main="LeadingEdge Heatmap GSE16008", col=bluered,

         row.clustering="supervised", col.clustering="supervised",

         row.dendrogram=TRUE, ColSideColors=gluc.class)


dev.off()

###Boxplot SmokingCessation

##########################

# gluc.class <- c(rep("1-Baseline",8), rep("2-four", 6), rep("3-eight", 7), rep("4-sixteen", 6), rep("5-twentyfour", 6))
# m2 <- lm(train.scores~as.numeric(Time1)) 
# pval<-anova(m2)[2,4]
# pdf(file="/protected/projects/pulmarray/U01_SmokingCessation/Paper/Figures/Figure3B Boxplot of leading edge genes GSE16008 SmokingCessation.pdf")

 
	# boxplot(as.numeric(-train.scores)~gluc.class, xaxt="n",yaxt="n",col=rep("cyan",length(vsamples)),main=pval,ylab='Metagene Score',xlab='Time in Weeks')

	# axis(1,at=c(1:5), labels=c("0","4","8","16","24"))
	# axis(2,at=NULL, labels=TRUE,tick=TRUE)

# dev.off()
