################################################Normalization and Entrez to gene symbol conversion##############################################################
#Affy																																						   #
               																																					   #
################################################################################################################################################################
	#R version 2.13.1
	library(affy)
    library(hugene10stv1hsentrezgcdf, lib.loc="/data/share/pulm/BrainArray/14.1.0")                                                        ##Calling cdf version 14.1.0
    setwd("/protected/projects/pulmarray/U01_SmokingCessation/Cell_Files_All/")
    exprs<-c()
    exprs<-justRMA(cdfname="HuGene10stv1_Hs_ENTREZG")                                          ##get RMA data
	
	
	write.exprs(exprs,"/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/SC_rmaexpfile_Entrez_14.1.0.txt")                                        ##Writes expression values to text file in working directory for affy package
	
################################################Linear mixed effects model, heatmaps, boxplots###############################################################
#																																							#
#																																							   #
#############################################################################################################################################################
####Final Program For Time taken as 0481624##############
setwd("/protected/projects/pulmarray/U01_SmokingCessation/")
Exp<-read.delim("Paper/Files/SC_rmaexpfile_Entrez_14.1.0.txt",check.names="FALSE")   				##Read the entrez id expression file
                                    				
gene.id <- Exp[,1]
Exp.x <- Exp[,2:58]
rownames(Exp.x) <- as.character(gene.id)
Exp <- Exp.x
                                    				            
DemoBmc <- read.table ("Paper/Files/BMC_SmokingCessation_Demographics.txt",sep="\t",header=TRUE) #Read Demographic file BMC,UM
DemoMn <- read.table ("Paper/Files/MN_SmokingCessation_Demographics.txt",sep="\t",header=TRUE) 
source("Paper/R_Programs_Paper/zscore_rows.R")
#new<-zscore.rows(as.matrix(log2(Exp))
data.pca <- prcomp(Exp),scale=TRUE,center=TRUE) #Doing PCA
#color <- factor(ifelse(phenoData$Annot_Cancer_Status=='Cancer','black','yellow'))
#shape <- factor(ifelse(phenoData$Annot_Cancer_Status=='Cancer',21,25))
plot(data.pca$rotation[,1],data.pca$rotation[,2], xlab="PC 1", ylab="PC 2")#main="", col=color,bg=color #plot

text(data.pca$rotation[,1] , data.pca$rotation[,2],pos=3,labels=as.character(colnames(Exp)))#put labels on the plot


patient.bmc <- c(107,116,206)		 									#Creating a vector of patients for BMC excluding outliers, outliers --> 162 (2 samples),213 (4 samples)
patient.mn <- c(205,208,227,243,265) 									#Creating a vector of patients for UM excluding outliers, outliers --> 219 (5 samples),258 (3 samples)

patient.ind.bmc=c()														# Rearranging the BMC demographic file excluding outliars
for(i in 1:length(patient.bmc)){
bmc<-grep(patient.bmc[i],DemoBmc[,1])
patient.ind.bmc <- append(patient.ind.bmc,bmc)}
DemoBmc1 <- DemoBmc[patient.ind.bmc,]									# Creating new demo file for BMC
DemoBMC.out<-DemoBmc[-patient.ind.bmc,]									# Creating demofile with outliers

patient.ind.mn=c()														# Rearranging the UM demographic file excluding outliars
for(i in 1:length(patient.mn)){
mn<-grep(patient.mn[i],DemoMn[,2])
patient.ind.mn <- append(patient.ind.mn,mn)}
DemoMn1 <- DemoMn[patient.ind.mn,]										# Creating new demo file for BMC
DemoMn.out<-DemoMn[-patient.ind.mn,]									# Creating demofile with outliers



cell.id<-DemoMn1[,1]													#Extracting UM cell name id's to match expression file
patient <- c(patient.bmc,cell.id)                   					#Creating a single vector of matching patient id's in Demographic and Exp file

patient.ind <- c()														#Extracting out indices for column names from Exp and matching against demographic files
for(i in 1:length(patient)){
ind<-grep(patient[i],colnames(Exp))
patient.ind <- append(patient.ind,ind)}
Exp2 <- Exp[,patient.ind]												#Creating new expression file with reduced number of columns / samples dim --> 19741 36



All_samples <- c(as.character(DemoBmc1[,1]),as.character(DemoMn1[,2]))  #Create a vector of all samples with time suffix

excl_samples <- c("N205_4"," N205_8  ","N227_16")                       #Samples as outliars

ex.ind <- match(excl_samples,All_samples,nomatch=0)                      #Get indices for EXp file
Mn.ind <-match(excl_samples,as.character(DemoMn1[,2]),nomatch=0)         #Get indices for UM file
bmc.ind<- match(excl_samples,as.character(DemoBmc1[,1]),nomatch=0)       #Get indices for BMC file
###Exclude samples
E<-Exp2[,-ex.ind]
B<-DemoBmc1
D <- DemoMn1[-Mn.ind,]
DemoMn2 <-D                                                             #Create new Demograohic file for UM

DemoBmc2<-B                                                             #Create new Demograohic file for BMC
Exp.ex <-E 																#Create new Exp file with further reduced columns/samples dim --> 19741 31

Time1 <- c(DemoBmc2[,3],DemoMn2[,4])                                     #Get vector of Time 
equal0<- which(Time1==0)
equal4 <- which(Time1==4)
equal8 <- which(Time1==8)
equal16 <- which(Time1==16)
equal24 <- which(Time1==24)
Times <- replace(Time1,equal0,0)
Times <- replace(Times,equal4,0)
Times <- replace(Times,equal8,1)
Times <- replace(Times,equal16,1)
Times <- replace(Times,equal24,1)

patient_id <- c(DemoBmc2[,5],DemoMn2[,3])							    #patient id
patient_id <-as.factor(patient_id)
RIN<-c(DemoBmc2[,13],DemoMn2[,9])									  # Get RIN scores
RIN<-as.numeric(RIN)


res <-c()
library(nlme)
for(i in 1:nrow(Exp.ex)){

gene.exp<-as.numeric(Exp.ex[i,])
try({

m1<-lme(gene.exp~as.factor(Times)+RIN,random=~1|patient_id,method="ML")
m2<-lme(gene.exp~RIN,random=~1|patient_id,method="ML")

#calculate and see the p-value for factor y
coef<-coef(m1)[1,2]
tstat<-summary(m1)$tTable[2,4]
slope <- summary(m1)$tTable[2,1]
pval<-anova(m1,m2)[2,9]	

cp<-cbind(as.character(gene.id[i]),coef,tstat,slope,pval)

print (i)
res<-rbind(res,cp)

},TRUE)
rownames(res)<-c()

}

rowname <- res[,1]
res1 <- res[,2:5]
rownames(res1) <- rowname
res.fdr<-p.adjust(as.numeric(res1[,4]),method="fdr")
res2<-cbind(res1,res.fdr)
write.table(res2,"SC_outputmodel_results.txt",sep="\t")
#res2<- read.table("SC_outputmodel_results.txt",sep="\t")
index <- match(rownames(res2),as.character(gene.id),nomatch=0)                         # Match gene id that passed the model
Exp.ex1 <- Exp.ex[index,]    
gs <- gene.id[index]                                                                   # Create matrice with matched gene id
rownames(Exp.ex1) <- gs

annot<-read.table("Paper/Files/ModelResults/062713_SC_annotationFile_Entrez_Gene_allgenes_cdf_14.1.0.txt",sep="\t")
na <- which(annot[,1]!="NA") #19741
annot1 <- annot[na,]
annot <-annot1 #substracting out dataframe of na 19197
list_entrez <- rownames(res2)
entrezid<-rownames(annot)
Entrezexp <- as.character(list_entrez)
entrezid <- as.character(entrezid)
exp_ind<-match(entrezid,Entrezexp,nomatch=0)
Entrex_ind<-match(Entrezexp,entrezid,nomatch=0)
gene.symbol<-as.character(annot[Entrex_ind,1])
res3 <- res2[exp_ind,]
Exp.ex2 <- Exp.ex1[exp_ind,]
rownames(Exp.ex2)<-c()
rownames(Exp.ex2) <- gene.symbol
res.new <-cbind(as.character(gene.symbol),res3)


fdr05<-which(res.new[,6]<0.05)                                                      #Select genes less then fdr 05
length(fdr05)   
out<- res.new[fdr05,]
Time.ord <- order(Time1)                                                             #Order time increasing
Time2 <- Time1[Time.ord]                                                             #Adjust Time according to order
NewExp <- Exp.ex2[fdr05,Time.ord]													#Order all required data in time order and significant genes 3314 33
patient_id1 <- patient_id[Time.ord]
RIN1 <- RIN[Time.ord]

down<-which(as.numeric(out[,2])<0) 																    #Getting indices for negative FC
up<-which(as.numeric(out[,2])>=0) 															    #Getting indices for positive FC
dfc.ind <-which(abs(as.numeric(out[down,2]))<8)  														    #Calculating and storing linear FC for up indices                                                         
ufc.ind<-which(as.numeric(out[up,2])>=8)  											   #Calculating and storing linear FC for down indices
  											   #Selecting up genes with FC >1.7
length(dfc.ind)                                                                    #
length(ufc.ind)                                                                    #
gene.symbol.gb1<-out[c(down[dfc.ind],up[ufc.ind]),]
#gene.symbol.gb1<-out
gene.symbol.gb <- gene.symbol.gb1[,1]								           # genes up and down at fdr and FC cutoff

j<-match(gene.symbol.gb,rownames(NewExp),nomatch=0)
NewExp1<-NewExp[j,]

write.table(gene.symbol.gb1,"SmokingCess_00111C_Genes_fdr05_2188D_7U_genesymbol_results.txt",sep="\t")
write.table(gene.symbol.gb1,"SmokingCess_00111C_Genes_fdr05_genesymbol_results.txt",sep="\t")
write.table(gene.symbol.gb1,"SmokingCess_00011C_Genes_fdr05_D_U_genesymbol_results.txt",sep="\t")
write.table(gene.symbol.gb1,"SmokingCess_00001C_Genes_fdr05_D_U_genesymbol_results.txt",sep="\t")
write.table(gene.symbol.gb1,"SmokingCess_Timenumeric_Genes_fdr05_D_U_genesymbol_results.txt",sep="\t")


####check for overlaps

##################################non linear modelsoverlap####################
file<-read.delim("Paper/Files/ModelResults/SC_LME_Time_Entrez_genesym_Model_Results.txt")
file1<-file[which(file[,5]<0.05),]
file2<-read.delim("SmokingCess_00111C_Genes_fdr05_2188D_7U_genesymbol_results.txt")

file3<-read.delim("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/ModelResults/SC_Genes_fdr05_FC1.7_96D_23U_genesymbol_results.txt")
file4<-read.delim("SmokingCess_00111C_Genes_fdr05_genesymbol_results.txt")

length(intersect(file3[,1],file2[,1])) #113 common at fdr and FC

length(intersect(file1[,1],file4[,1])) #3206 
#############################calculating FC using gene expression data
NewExp1<-as.data.frame(NewExp)
mean0<-mean(as.numeric(NewExp1[1,which(Times==0)]))
mean1<-mean(as.numeric(NewExp1[1,which(Times==1)]))
Exp0<-NewExp1[,which(Times==0)]
Exp1<-NewExp1[,which(Times==1)]
mean0<-apply(as.matrix(Exp0),1,mean)
mean1<-apply(as.matrix(Exp1),1,mean)
log2FC<-mean0-mean1
FC<-2^log2FC
dn<-which(log2FC<0.766)#2084
up<-which(log2FC>=0.766)#4289
slopes.FC.up <- FC[up]  														    #Calculating and storing linear FC for up indices                                                         
slopes.FC.down<-FC[dn]												   #Calculating and storing linear FC for down indices

dfc.ind <- which(slopes.FC.down < 1.7)                                               #Selecting down genes with FC >1.7																			
ufc.ind <- which(slopes.FC.up>=1.7)    											   #Selecting up genes with FC >1.7
length(dfc.ind)                                                                    #
length(ufc.ind)                                                                    #
gene.symbol.gb1<-out[c(down[dfc.ind],up[ufc.ind]),]
gene.symbol.gb <- gene.symbol.gb1[,1]	
######convert mouse genesets to human genesets
setwd("/protected/projects/pulmarray/U01_SmokingCessation/")
ss<-read.delim("Mice_CS_vs_SS.txt")
cs<-read.delim("Mice_CS_vs_AC.txt")
convertEnsHumanGSMouse <- function(x){
 
require("biomaRt")
human = useMart(host="feb2014.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
mouse = useMart(host="feb2014.archive.ensembl.org",biomart="ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = "mgi_symbol", filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL =human , uniqueRows=T)
humanx <- unique(genesV2[, 2])
 
# Print the first 6 genes found to the screen
print(head(humanx))
#return(humanx)
return(genesV2)
}
ssH<-convertEnsHumanGSMouse(ss[,1])

ssH1<-ssH[match(ss[,1],ssH[,1],nomatch=0),]
ss1<-ss[match(ssH1[,1],ss[,1],nomatch=0),]
newSSH<-cbind(ssH1[,2],ss1[,2])
write.table(newSSH,"HG_Mice_CS_vs_SS.txt",sep="\t")

csH<-convertEnsHumanGSMouse(cs[,1])
csH1<-csH[match(cs[,1],csH[,1],nomatch=0),]
cs1<-cs[match(csH1[,1],cs[,1],nomatch=0),]
newCSH<-cbind(csH1[,2],cs1[,2])
write.table(newCSH,"HG_Mice_CS_vs_AC.txt",sep="\t")
#########residuals
Matrix<-as.matrix(NewExp1)														  #Adjust residuals for patient effect and RIN

name<-colnames(Matrix)
model <- model.matrix(~1 + patient_id1+RIN1)									  #Get residuals after adjusting for patient and RIN

fileGE<-as.matrix(Matrix) %*% (diag(dim(Matrix)[2]) - model %*% solve(t(model) %*% model) %*% t(model)) #Get residuals after adjusting for patient and RIN

colnames(fileGE) <- colnames(Matrix)


colnames(fileGE)<-c()
names<-c(as.character(DemoBmc2[,1]),as.character(DemoMn2[,2]))
colnames(fileGE)<-names[Time.ord]                                                         #Change colnames to simpler names


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
P1<- fileGE[,1:8]
P1<-cbind(P1,c(rep(0,length(index))))
P2<- fileGE[,9:14]
P2<-cbind(P2,c(rep(0,length(index))))
P3<- fileGE[,15:21]
P3<-cbind(P3,c(rep(0,length(index))))
P4<- fileGE[,22:27]
P4<-cbind(P4,c(rep(0,length(index))))
P5<- fileGE[,28:33]

E<-cbind(P1,P2,P3,P4,P5)
source("Paper/R_Programs_Paper/zscore_rows.R")
new<-zscore.rows(E)
#################################################Heatmap###################################################
vsamples<-c(8,1,6,1,7,1,6,1,6)
library(heatmap3)                                            #Load graphics/plotting packages that has red/blue colors for the heat map

library(gplots)                                              #create a variable with colors that will be used in heat map

bluered<-bluered(256)
##Creating a vector that contains colorsrepresenting the class of each sample
Tclasslables<-c(rep("darkblue",vsamples[1]),rep("white",vsamples[2]),rep("blue",vsamples[3]),rep("white",vsamples[4]),
		rep("darkgreen",vsamples[5]),rep("white",vsamples[6]),rep("pale green",vsamples[7]),rep("white",vsamples[8]),rep("yellow",vsamples[9]))
		
		

pdf(file="Paper/Figures/SuppleFigure_SC_LME_8patients_33samples_114 genes FDR05+FC 1.7_00111_Heatmap.pdf")	


heatmap3(E, main="Heatmap of 114 genes adjusted for patient and rin effect", col=bluered,

         row.clustering="unsupervised", col.clustering="supervised",

         row.dendrogram=TRUE, ColSideColors=Tclasslables)


dev.off()	                                                                       #close output file
