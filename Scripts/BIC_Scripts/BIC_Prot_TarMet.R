#!/usr/bin/Rscript
index <- as.integer(Sys.getenv("PBS_ARRAYID"))
library("bnlearn")

#	Load Data
relat_ProtTarMet <- read.table("/users/home/dafmic/CrossSectional//RELATIONS/results_TOTALS/ProtTarMet_Assoc_Signif.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
relat_ProtTarMet$MetaboliteID <- trimws(relat_ProtTarMet$MetaboliteID)
covar_Protf <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Proteomics/Proteins_COVARIATES.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
raw_proteins_Imp <-  read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Proteomics/Proteins_processed.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_proteins_TechBiol <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/Prot_Residuals_TechBiol_withAgeEffects.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_proteins_TechBiol$ProteinID <-trimws(resid_proteins_TechBiol$ProteinID)
resid_tarMet_TechBiol <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/TarMet_TechBiol_withAgeEffects.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_tarMet_TechBiol_sub3023<- data.frame(resid_tarMet_TechBiol[,1],resid_tarMet_TechBiol[,c(which(colnames(resid_tarMet_TechBiol) %in% colnames(raw_proteins_Imp)))])
colnames(resid_tarMet_TechBiol_sub3023)[2:3024] <- colnames(raw_proteins_Imp)
colnames(resid_tarMet_TechBiol_sub3023)[1] <- "Metabolite"
resid_tarMet_TechBiol_sub3023$Metabolite <- trimws(resid_tarMet_TechBiol_sub3023$Metabolite)
#*******************************************************************************************		 PROT -- TARMET	******************************************************************************************************
modelINDEP <- empty.graph(c("Age", "Protein", "TarMet"))
modelTARMET <- empty.graph(c("Age", "Protein", "TarMet"))
modelPROTEIN <- empty.graph(c("Age", "Protein", "TarMet"))
modelstring(modelINDEP) <- "[Protein][TarMet][Age][Protein|Age][TarMet|Age]" # Independet model
modelstring(modelTARMET) <- "[Protein][TarMet][Age][Protein|Age][TarMet|Protein]" # Protein via TarMet
modelstring(modelPROTEIN) <- "[Protein][TarMet][Age][TarMet|Age][Protein|TarMet]" # TarMet via Protein
#	19,853
BIC_Prot_TarMet <- matrix(0,ncol=5,nrow=1986)
colnames(BIC_Prot_TarMet) <- c("ProteinID","MetaboliteID","modelPROTEIN", "modelTARMET","modelINDEP")
relat_ProtTarMet2 <- relat_ProtTarMet[(index*1986-1985):min(index*1986,nrow(relat_ProtTarMet)),]

for(i in 1:nrow(relat_ProtTarMet2)){
	ProteinID <- as.character(relat_ProtTarMet2[i,1])
	MetaboliteID <- as.character(relat_ProtTarMet2[i,3])
	age <- as.numeric(covar_Protf$Age)
	index.p <- which(as.data.frame(resid_proteins_TechBiol)$ProteinID == ProteinID)
	index.m <- which(as.data.frame(resid_tarMet_TechBiol_sub3023)$Metabolite == MetaboliteID)
	p <- as.numeric(resid_proteins_TechBiol[index.p,3:3025])
	m <- as.numeric(resid_tarMet_TechBiol_sub3023[index.m,2:3024])
	one_trio <- data.frame(p,m,age)
	colnames(one_trio)[1] <- 'Protein'
	colnames(one_trio)[2] <- 'TarMet'
	colnames(one_trio)[3] <- 'Age'
	s1 <- score(modelPROTEIN, data = one_trio, type ='bic-g')
    s2 <- score(modelTARMET, data = one_trio, type ='bic-g')
    s3 <- score(modelINDEP, data = one_trio, type ='bic-g')
    BIC_Prot_TarMet[i,] <- c(ProteinID,MetaboliteID,s1,s2,s3)
	print(i)
}

BIC_Prot_TarMet <- data.frame(BIC_Prot_TarMet)
for(i in 1:nrow(BIC_Prot_TarMet)){
	s1 <- as.numeric(BIC_Prot_TarMet[i,"modelPROTEIN"])
	s2 <- as.numeric(BIC_Prot_TarMet[i,"modelTARMET"])
	s3 <- as.numeric(BIC_Prot_TarMet[i,"modelINDEP"])
	if(abs(s1-s2) > (s1-s3)){ # s2 < s3
		if(s3 < s1){	# s2 < s3 < s1
			BIC_Prot_TarMet$TOP[i] <- "modelTARMET" 
			BIC_Prot_TarMet$DIFF[i] <- abs(s2-s3)}	
		else if(s1 < s2){	# s1 < s2 < s3
			BIC_Prot_TarMet$TOP[i] <- "modelPROTEIN" 
			BIC_Prot_TarMet$DIFF[i] <- abs(s1-s2)}
		else{ # s2 < s1 < s3
			BIC_Prot_TarMet$TOP[i] <- "modelTARMET" 
			BIC_Prot_TarMet$DIFF[i] <- abs(s2-s1)}}
	else{									# s3 < s2
		if(s2 < s1){	# s3 < s2 < s1
			BIC_Prot_TarMet$TOP[i] <- "modelINDEP" 
			BIC_Prot_TarMet$DIFF[i] <- abs(s3-s2)}
		else if(s1 < s3){ 	# s1 < s3 < s2
			BIC_Prot_TarMet$TOP[i] <- "modelPROTEIN" 
			BIC_Prot_TarMet$DIFF[i] <- abs(s1-s3)}
		else{		# s3 < s1 < s2
			BIC_Prot_TarMet$TOP[i] <- "modelINDEP" 
			BIC_Prot_TarMet$DIFF[i] <- abs(s3-s1)}}
}
# The difference between the BICs for two models (∆BIC) can be computed,
# where guidelines indicate that a difference greater than 10 indicates very strong evidence;
# 6–10 indicates strong evidence;
# 2–6 indicates positive evidence;
# and 0–2 indicates weak evidence for the more complex model (Raftery, 1995)
for(i in 1:nrow(BIC_Prot_TarMet)){
	if(BIC_Prot_TarMet$DIFF[i] > 10){
		BIC_Prot_TarMet$SIGNIF[i]  <- "***" }
	else if(BIC_Prot_TarMet$DIFF[i]>6 && BIC_Prot_TarMet$DIFF[i]<10){
	BIC_Prot_TarMet$SIGNIF[i]  <- "**" }
	else if(BIC_Prot_TarMet$DIFF[i]>2 && BIC_Prot_TarMet$DIFF[i] <6){
	BIC_Prot_TarMet$SIGNIF[i]  <- "*" }
	else {BIC_Prot_TarMet$SIGNIF[i]  <- "-" }	
}
write.table(BIC_Prot_TarMet, file = paste("/users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Prot_TarMet", index, sep="_"), quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
#colnames(BIC_Prot_TarMet) <- c("ProteinID","MetaboliteID","modelPROTEIN", "modelTARMET","modelINDEP","TOP","DIFF","SIGNIF")