
##Figure 6 DAVID catagories boxplots
###################################################SmokingCessation data################################################
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


Time <- c(DemoBmc2[,3],DemoMn2[,4])

time.ind <-order(Time)
Time.ord <- Time[time.ind]
RIN.ord <- RIN[time.ind]
patient_id.ord <- patient_id[time.ind]
NewExp<-Exp[,time.ind]
source("/protected/projects/pulmarray/U01_SmokingCessation/Paper/R_Programs_Paper/zscore_rows.R")
Exp2<-zscore.rows(NewExp)
		
			

###############################File with different functional categories and genes########################

func.file <- read.delim("David_functional_categories_list.txt")
func.xeno<-func.file[1:5,1] # Metabolism of xeno
func.apop<- func.file[1:9,2] #anti-apoptosis
func.hom <- func.file[1:23,3]#Homeostasis
func.wnd <- func.file[1:10,4]#response to wounding

#########match with expression file
xeno.ind<- match(as.character(func.xeno),genes,nomatch=0)
Exp.xeno <- Exp2[xeno.ind,]

xeno.pca<-prcomp(Exp.xeno,center=F)
pc.xeno<-xeno.pca$rotation[,1]

apop.ind<- match(as.character(func.apop),genes,nomatch=0)
Exp.apop <- Exp2[apop.ind,]
apop.pca<-prcomp(Exp.apop,center=F)
pc.apop<-apop.pca$rotation[,1]

hom.ind<- match(as.character(func.hom),genes,nomatch=0)
Exp.hom <- Exp2[hom.ind,]
hom.pca<-prcomp(Exp.hom,center=F)
pc.hom<-hom.pca$rotation[,1]


wnd.ind<- match(as.character(func.wnd),genes,nomatch=0)
Exp.wnd <- Exp2[wnd.ind,]
wnd.pca<-prcomp(Exp.wnd,center=F)
pc.wnd<-wnd.pca$rotation[,1]

All.pc<- cbind(pc.xeno,pc.apop,pc.hom,pc.wnd)
All.pc.t <- t(All.pc)

gluc.class <- c(rep("1-Baseline",8), rep("2-four", 6), rep("3-eight", 7), rep("4-sixteen", 6), rep("5-twentyfour", 6))
library(nlme)

list.paths<- c("Metabolism of Xeno","Apoptosis","Homeostasis","ResponseToWounding")
pdf(file="/protected/projects/pulmarray/U01_SmokingCessation/Paper/Figures/Figure6B David functional categories1.pdf",title="BoxPlots of functional categories")
par(mar=rep(2,4))
layout(matrix(1:4,nrow=2,ncol=2,byrow=TRUE))

for (i in 1:nrow(All.pc.t)){
pc1<- All.pc.t[i,]
pval<-c()
m1<-lme(pc1~as.numeric(Time.ord)+as.numeric(RIN.ord),random=~1|patient_id.ord,method="ML")
m2<-lme(pc1~as.numeric(RIN.ord),random=~1|patient_id.ord,method="ML")

#calculate and see the p-value 

	pval<-anova(m1,m2)[2,9]	
	print(pval)
	
	boxplot(as.numeric(-pc1)~gluc.class, xaxt="n",yaxt="n",ylab='PC1 GeneExp',xlab='Time in Months')

	axis(1,at=c(1:5), labels=c("0","4","8","16","24"))
	axis(2,at=NULL, labels=TRUE,tick=TRUE)
	



}
dev.off()

pdf(file="/protected/projects/pulmarray/U01_SmokingCessation/Paper/Figures/test.pdf")
	boxplot(as.numeric(pc.apop)~gluc.class, xaxt="n",yaxt="n",col=rep("light green",length(vsamples)),main=c(list.paths[2],pval),ylab='PC1 GeneExp',xlab='Time in Months')

dev.off()