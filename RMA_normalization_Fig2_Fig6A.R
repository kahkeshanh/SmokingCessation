################################################Normalization and Entrez to gene symbol conversion##############################################################
#Affy																																						   #
               																																					   #
################################################################################################################################################################
	library(affy);
    library(affyio);
    library(affyPLM);
	
	setwd("/protected/projects/pulmarray/U01_SmokingCessation/Cell_Files_All/")
	output.path="/protected/projects/pulmarray/U01_SmokingCessation"
	output.prefix <- "SmokignCessation_2019" # e.g., "Nasaldiseases_180502"

	x<-read.delim("list.txt")

cat("Computing RLE and NUSE metrics.\n");
abatch <- ReadAffy()
pset <- fitPLM(abatch)
QC <- list();
QC$RLE <- RLE(pset, type="values");
QC$NUSE <- NUSE(pset, type="values");

QC.medians <- list();
for (metric in names(QC)) {
    QC.medians[[metric]] <- apply(QC[[metric]], 2, median);
}

# Save matrices of QC metrics as RDS files
for (metric in names(QC)) {
    output.filename <- file.path(output.path, sprintf("%s_%s.rds", output.prefix, metric));
    saveRDS(QC[[metric]], output.filename);
}
QC.cutoffs <- list(RLE = 0.2, NUSE = 1.05);	

###
match(colnames(QC$RLE,)
####

	output.filename <- file.path(output.path, sprintf("%s_RLE_NUSE.pdf", output.prefix));
    #cat(sprintf("Drawing RLE and NUSE boxplots to '%s'.\n", output.filename));
    pdf("/protected/projects/pulmarray/U01_SmokingCessation/RLENUSE.pdf", width=11, height=8.5);
    boxplot(
        QC$RLE,
        main="RLE (Relative Log Expression)\nShould be centered on 0 (blue line)\nRed sample = out of bounds (dashed red line)",
        names=colnames(QC$RLE), las=2,
        border=c("black","red")[(QC.medians$RLE > QC.cutoffs$RLE)+1]    
    );
    lines(x=c(0,ncol(QC$RLE)+1), y=rep(0,2), col="blue", lty=2);
    lines(x=c(0,ncol(QC$RLE)+1), y=rep(QC.cutoffs$RLE,2), col="red", lty=2);
    boxplot(
        QC$NUSE,
        main="NUSE (Normalized Unscaled Standard Error)\nShould be centered on 1 (blue line)\nRed sample = out of bounds (dashed red line)",
        names=colnames(QC$NUSE), las=2
        #border=c("black","red")[(QC.medians$NUSE > QC.cutoffs$NUSE)+1]  
    );
    lines(x=c(0,ncol(QC$NUSE)+1), y=rep(1,2), col="blue", lty=2);
    lines(x=c(0,ncol(QC$NUSE)+1), y=rep(QC.cutoffs$NUSE,2), col="red", lty=2);
    dev.off();

	#R version 2.13.1
	library(affy)
    library(hugene10stv1hsentrezgcdf, lib.loc="/data/share/pulm/BrainArray/14.1.0")                                                        ##Calling cdf version 14.1.0
    
    exprs<-c()
    exprs<-justRMA(cdfname="HuGene10stv1_Hs_ENTREZG")                                          ##get RMA data
	
	
	write.exprs(exprs,"/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/SC_rmaexpfile_Entrez_14.1.0.txt")                                        ##Writes expression values to text file in working directory for affy package
	
################################################Linear mixed effects model, heatmaps, boxplots###############################################################
#																																							#
#																																							   #
#############################################################################################################################################################
####Final Program For Time taken as 0481624##############
setwd("/protected/projects/pulmarray/U01_SmokingCessation/")
Exp<-read.delim("Paper/Files/SC_rmaexpfile_Entrez_14.1.0.txt")   				##Read the entrez id expression file
                                 				
gene.id <- Exp[,1]
Exp.x <- Exp[,2:58]
rownames(Exp.x) <- as.character(gene.id)
Exp <- Exp.x
                                    				            
DemoBmc <- read.table ("Paper/Files/BMC_SmokingCessation_Demographics.txt",sep="\t",header=TRUE) #Read Demographic file BMC,UM
DemoMn <- read.table ("Paper/Files/MN_SmokingCessation_Demographics.txt",sep="\t",header=TRUE) 


		
	


patient.bmc <- c(107,116,206)		 									#Creating a vector of patients for BMC excluding outliers, outliers --> 162 (2 samples),213 (4 samples)
patient.mn <- c(219,208,227,243,265) 									#Creating a vector of patients for UM excluding outliers, outliers --> 219 (5 samples),258 (3 samples)

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

Time <- c(DemoBmc2[,3],DemoMn2[,4])                                     #Get vector of Time 
Time <- Time+1
Time1 <- log2(Time)														#log transform time
Time1 <- as.numeric(Time1)
patient_id <- c(DemoBmc2[,5],DemoMn2[,3])							    #patient id
patient_id <-as.factor(patient_id)
RIN<-c(DemoBmc2[,13],DemoMn2[,9])									  # Get RIN scores
RIN<-as.numeric(RIN)


res <-c()
library(nlme)
for(i in 1:nrow(Exp.ex)){

gene.exp<-as.numeric(Exp.ex[i,])
try({

m1<-lme(gene.exp~Time1+RIN,random=~1|patient_id,method="ML")
m2<-lme(gene.exp~RIN,random=~1|patient_id,method="ML")

#calculate and see the p-value for factor y
coeff<-summary(m1)$tTable[2,4]
slope <- summary(m1)$tTable[2,1]
pval<-anova(m1,m2)[2,9]	

cp<-cbind(as.character(gene.id[i]),coeff,slope,pval)

print (i)
res<-rbind(res,cp)

},TRUE)
rownames(res)<-c()

}


rowname <- res[,1]
res1 <- res[,2:4]
rownames(res1) <- rowname
res.fdr<-p.adjust(as.numeric(res1[,3]),method="fdr")
res2<-cbind(res1,res.fdr)
#write.table(res2,"SC_outputmodel_results.txt",sep="\t")
#res2<- read.table("SC_outputmodel_results.txt",sep="\t")
index <- match(rownames(res2),as.character(gene.id),nomatch=0)                         # Match gene id that passed the model
Exp.ex1 <- Exp.ex[index,]    
gs <- gene.id[index]                                                                   # Create matrice with matched gene id
rownames(Exp.ex1) <- gs
####Convert entrez-->gene symbols#############################################################################################
library(hugene10stv1hsentrezg.db , lib.loc="/data/share/pulm/BrainArray/14.1.0")
#Call for annotation scriptfrom other directory
source("/protected/projects/pulmarray/U01_SmokingCessation/Paper/R_Programs_Paper/probeset.annotation.R")
#Getting the list of Entrez gene id's
#Read in table with Entrez gene id's
list_entrez <- rownames(res2)
#call to function
annot <- probeset.annotation(list_entrez, platform="hugene10stv1", mapping="hsentrezg", fields=c("SYMBOL", "GENENAME"))
####Look for no matches and remove
na <- which(annot[,1]!="NA") #19741
annot1 <- annot[na,]
annot <-annot1 #substracting out dataframe of na 19197
#write.table(annot,"Paper/Files/ModelResults/062713_SC_annotationFile_Entrez_Gene_allgenes_cdf_14.1.0.txt",sep="\t")
annot<-read.delim("Paper/Files/ModelResults/062713_SC_annotationFile_Entrez_Gene_allgenes_cdf_14.1.0.txt")

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

#write.table(res.new,"Paper/Files/ModelResults/SC_LME_Time_Entrez_genesym_Model_Results.txt",sep="\t")	                                   #writing the result from the model
#write.table(Exp.ex2,"/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/ModelResults/SC_ExpressionData_All_genes_genesymbol.txt",sep="\t")
Exp.ex2 <- read.table("/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/ModelResults/SC_ExpressionData_All_genes_genesymbol.txt",sep="\t")
res.new<- read.table("Paper/Files/ModelResults/SC_LME_Time_Entrez_genesym_Model_Results.txt",sep="\t")
###########################################################################################################################

fdr05<-which(res.new[,5]<0.05)                                                      #Select genes less then fdr 05
length(fdr05)   #3452
res05<-res.new[fdr05,]

write.table(res05,"/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/ModelResults/SC_Genes_fdr05_genesymbol_Modelresults.txt",sep="\t")

out<- res.new[fdr05,]
Time.ord <- order(Time)                                                             #Order time increasing
Time2 <- Time[Time.ord]                                                             #Adjust Time according to order
NewExp <- Exp.ex2[fdr05,Time.ord]													#Order all required data in time order and significant genes 3314 33
patient_id1 <- patient_id[Time.ord]
RIN1 <- RIN[Time.ord]

FC <- as.numeric(out[,3])*4.643856                                                  #Calculate fold change by slope
out1<-cbind(out,FC)
down<-which(FC<0) 																    #Getting indices for negative FC
up<-which(FC>=0)  																    #Getting indices for positive FC
slopes.FC.up <- 2^FC[up]  														    #Calculating and storing linear FC for up indices                                                         
slopes.FC.down<-(2^abs(FC[down]))												   #Calculating and storing linear FC for down indices

dfc.ind <- which(slopes.FC.down > 1.7)                                               #Selecting down genes with FC >1.7																			
ufc.ind <- which(slopes.FC.up > 1.7)    											   #Selecting up genes with FC >1.7
length(dfc.ind)                                                                    #
length(ufc.ind)                                                                    #
gene.symbol.gb1<-out[c(down[dfc.ind],up[ufc.ind]),]
gene.symbol.gb <- gene.symbol.gb1[,1]								           # genes up and down at fdr and FC cutoff

j<-match(gene.symbol.gb,rownames(NewExp),nomatch=0)
NewExp1<-NewExp[j,]

#write.table(gene.symbol.gb1,"/protected/projects/pulmarray/U01_SmokingCessation/Paper/Files/ModelResults/SC_Genes_fdr05_FC1.7_96D_23U_genesymbol_results.txt",sep="\t")
##############calculating residuals
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
source("zscore_rows.R")
new<-zscore.rows(E)
#################################################Heatmap###################################################
vsamples<-c(8,1,6,1,7,1,6,1,6)
library(heatmap3)                                            #Load graphics/plotting packages that has red/blue colors for the heat map

library(gplots)                                              #create a variable with colors that will be used in heat map

bluered<-bluered(256)
##Creating a vector that contains colorsrepresenting the class of each sample
Tclasslables<-c(rep("darkblue",vsamples[1]),rep("white",vsamples[2]),rep("blue",vsamples[3]),rep("white",vsamples[4]),
		rep("darkgreen",vsamples[5]),rep("white",vsamples[6]),rep("pale green",vsamples[7]),rep("white",vsamples[8]),rep("yellow",vsamples[9]))
		
		

pdf(file="Paper/Figures/Figure2 SC_LME_8patients_33samples_119 genes FDR05+FC 1.7_Heatmap.pdf")	


heatmap3(E, main="Heatmap of 119 genes adjusted for patient and rin effect", col=bluered,

         row.clustering="unsupervised", col.clustering="supervised",

         row.dendrogram=TRUE, ColSideColors=Tclasslables)


dev.off()	                                                                       #close output file


#################################Box plots of significant genes associated with smoking cessation################
source("Paper/R_Programs_Paper/zscore_rows.R")
new<-zscore.rows(NewExp1)
data.pca<-prcomp(new,center=F)
pc1<-data.pca$rotation[,1]
#################make labels for boxplot#########

gluc.class <- c(rep("1-Baseline",8), rep("2-four", 6), rep("3-eight", 7), rep("4-sixteen", 6), rep("5-twentyfour", 6))

pdf(file="Paper/Figures/Figure 6A GenesKinetics_119_genes_boxplots.pdf")

#Rapid kinetics, light green, slow kinetics, cyan 
	boxplot(as.numeric(-pc1)~gluc.class, xaxt="n",yaxt="n",col=rep("light green",length(vsamples)),main="Genes Kinetics",ylab='Metagene Score',xlab='Time in Weeks')

	axis(1,at=c(1:5), labels=c("0","4","8","16","24"))
	axis(2,at=NULL, labels=TRUE,tick=TRUE)

dev.off()

