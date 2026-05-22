#!/usr/bin/Rscript
index <- as.integer(Sys.getenv("PBS_ARRAYID"))
library("bnlearn")

#	Load Data
relat_ProtProt <- read.table("/users/home/dafmic/CrossSectional//RELATIONS/results_TOTALS/ProtProt_Assoc_Signif.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
relat_ProtProt$ProteinID2 <- trimws(relat_ProtProt$ProteinID2)
covar_Protf <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Proteomics/Proteins_COVARIATES.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_proteins_TechBiol <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/Prot_Residuals_TechBiol_withAgeEffects.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_proteins_TechBiol$ProteinID <- trimws(resid_proteins_TechBiol$ProteinID)
#*******************************************************************************************		PROT --- PROT		******************************************************************************************************
modelINDEP <- empty.graph(c("Age", "Prot1", "Prot2"))
modelPROT1 <- empty.graph(c("Age", "Prot1", "Prot2"))
modelPROT2 <- empty.graph(c("Age", "Prot1", "Prot2"))
modelstring(modelINDEP) <- "[Prot1][Prot2][Age][Prot1|Age][Prot2|Age]" # Independet model
modelstring(modelPROT1) <- "[Prot1][Prot2][Age][Prot1|Age][Prot2|Prot1]" # Prot1 via Prot2
modelstring(modelPROT2) <- "[Prot1][Prot2][Age][Prot2|Age][Prot1|Prot2]" # Prot2 via Prot1
#	131,052
BIC_Prot_Prot <- matrix(0,ncol=5,nrow=1311)
colnames(BIC_Prot_Prot) <- c("ProteinID1","ProteinID2","modelPROT1", "modelPROT2","modelINDEP")
relat_ProtProt2 <- relat_ProtProt[(index*1311-1310):min(index*1311,nrow(relat_ProtProt)),]

for(i in 1:nrow(relat_ProtProt2)){
	ProtID1 <- as.character(relat_ProtProt2[i,1])
	ProtID2 <- as.character(relat_ProtProt2[i,3])
	age <- as.numeric(covar_Protf$Age)
	index.g1 <- which(resid_proteins_TechBiol$ProteinID == ProtID1)
	index.g2 <- which(resid_proteins_TechBiol$ProteinID == ProtID2)
	g1 <- as.numeric(resid_proteins_TechBiol[index.g1,3:3025])
	g2 <- as.numeric(resid_proteins_TechBiol[index.g2,3:3025])
	one_trio <- data.frame(g1,g2,age)
	colnames(one_trio)[1] <- 'Prot1'
	colnames(one_trio)[2] <- 'Prot2'
	colnames(one_trio)[3] <- 'Age'
	s1 <- score(modelPROT1, data = one_trio, type ='bic-g')
    s2 <- score(modelPROT2, data = one_trio, type ='bic-g')
    s3 <- score(modelINDEP, data = one_trio, type ='bic-g')
    BIC_Prot_Prot[i,] <- c(ProtID1,ProtID2,s1,s2,s3)
}
###	FOR t=21 AND t=75
BIC_Prot_Prot <- data.frame(BIC_Prot_Prot)
for(i in 1:nrow(BIC_Prot_Prot)){
	s1 <- as.numeric(BIC_Prot_Prot[i,"modelPROT1"])
	s2 <- as.numeric(BIC_Prot_Prot[i,"modelPROT2"])
	s3 <- as.numeric(BIC_Prot_Prot[i,"modelINDEP"])
	if(BIC_Prot_Prot[i,3] == BIC_Prot_Prot[i,4]){
		BIC_Prot_Prot$TOP[i] <- "1-1"
		BIC_Prot_Prot$DIFF[i] <- 1}
	else if(abs(s1-s2) > (s1-s3)){ # s2 < s3
		if(s3 < s1){	# s2 < s3 < s1
			BIC_Prot_Prot$TOP[i] <- "modelPROT2" 
			BIC_Prot_Prot$DIFF[i] <- abs(s2-s3)}	
		else if(s1 < s2){	# s1 < s2 < s3
			BIC_Prot_Prot$TOP[i] <- "modelPROT1" 
			BIC_Prot_Prot$DIFF[i] <- abs(s1-s2)}
		else{ # s2 < s1 < s3
			BIC_Prot_Prot$TOP[i] <- "modelPROT2" 
			BIC_Prot_Prot$DIFF[i] <- abs(s2-s1)}}
	else{									# s3 < s2
		if(s2 < s1){	# s3 < s2 < s1
			BIC_Prot_Prot$TOP[i] <- "modelINDEP" 
			BIC_Prot_Prot$DIFF[i] <- abs(s3-s2)}
		else if(s1 < s3){ 	# s1 < s3 < s2
			BIC_Prot_Prot$TOP[i] <- "modelPROT1" 
			BIC_Prot_Prot$DIFF[i] <- abs(s1-s3)}
		else{		# s3 < s1 < s2
			BIC_Prot_Prot$TOP[i] <- "modelINDEP" 
			BIC_Prot_Prot$DIFF[i] <- abs(s3-s1)}}
}
# The difference between the BICs for two models (∆BIC) can be computed,
# where guidelines indicate that a difference greater than 10 indicates very strong evidence;
# 6–10 indicates strong evidence;
# 2–6 indicates positive evidence;
# and 0–2 indicates weak evidence for the more complex model (Raftery, 1995)
for(i in 1:nrow(BIC_Prot_Prot)){
	if(BIC_Prot_Prot$DIFF[i] > 10){
		BIC_Prot_Prot$SIGNIF[i]  <- "***" }
	else if(BIC_Prot_Prot$DIFF[i]>6 && BIC_Prot_Prot$DIFF[i]<10){
	BIC_Prot_Prot$SIGNIF[i]  <- "**" }
	else if(BIC_Prot_Prot$DIFF[i]>2 && BIC_Prot_Prot$DIFF[i] <6){
	BIC_Prot_Prot$SIGNIF[i]  <- "*" }
	else {BIC_Prot_Prot$SIGNIF[i]  <- "-" }	
}
write.table(BIC_Prot_Prot, file = paste("/users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Prot_Prot", index, sep="_"), quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
#colnames(BIC_Prot_Prot) <- c("ProteinID1","ProteinID2","modelPROT1", "modelPROT2","modelINDEP","TOP","DIFF","SIGNIF")