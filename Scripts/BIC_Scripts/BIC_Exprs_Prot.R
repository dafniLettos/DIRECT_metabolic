#!/usr/bin/Rscript
index <- as.integer(Sys.getenv("PBS_ARRAYID"))
library("bnlearn")

#	Load Data
relat_ExprsProt <- read.table("/users/home/dafmic/CrossSectional//RELATIONS/results_TOTALS/ExprsProt_Assoc_Signif.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
relat_ExprsProt$GeneID <- trimws(relat_ExprsProt$GeneID)
relat_ExprsProt$ProteinID <- trimws(relat_ExprsProt$ProteinID)
covar_Protf <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Proteomics/Proteins_COVARIATES.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
raw_proteins_Imp <-  read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Proteomics/Proteins_processed.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_genes_TechBiol <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/Exprs_Residuals_TechBiol_withAgeEffects.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_genes_TechBiol$GeneID <- trimws(resid_genes_TechBiol$GeneID)
resid_proteins_TechBiol <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/Prot_Residuals_TechBiol_withAgeEffects.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_proteins_TechBiol$ProteinID <- trimws(resid_proteins_TechBiol$ProteinID)
resid_genes_TechBiol_sub3023 <- data.frame(resid_genes_TechBiol[,1:2],resid_genes_TechBiol[,c(which(colnames(resid_genes_TechBiol) %in% colnames(raw_proteins_Imp)))])
colnames(resid_genes_TechBiol_sub3023)[3:3025] <- colnames(raw_proteins_Imp)
#*******************************************************************************************		EXPRS --- PROT		******************************************************************************************************
modelINDEP <- empty.graph(c("Age", "Exprs", "Protein"))
modelPROTEIN <- empty.graph(c("Age", "Exprs", "Protein"))
modelGENE <- empty.graph(c("Age", "Exprs", "Protein"))
modelstring(modelINDEP) <- "[Exprs][Protein][Age][Exprs|Age][Protein|Age]" # Independet model
modelstring(modelPROTEIN) <- "[Exprs][Protein][Age][Exprs|Age][Protein|Exprs]" # Gene via Protein
modelstring(modelGENE) <- "[Exprs][Protein][Age][Protein|Age][Exprs|Protein]" # Protein via Gene
#	1,036,751
BIC_Exprs_Prot <- matrix(0,ncol=5,nrow=1037)
colnames(BIC_Exprs_Prot) <- c("GeneID","ProteinID","modelGENE", "modelPROTEIN","modelINDEP")
relat_ExprsProt2 <- relat_ExprsProt[(index*1037-1036):min(index*1037,nrow(relat_ExprsProt)),]

for(i in 1:nrow(relat_ExprsProt2)){
	GeneID <- as.character(relat_ExprsProt2[i,1])
	ProteinID <- as.character(relat_ExprsProt2[i,3])
	age <- as.numeric(covar_Protf$Age)
	index.g <- which(resid_genes_TechBiol_sub3023$GeneID == GeneID)
	index.p <- which(as.data.frame(resid_proteins_TechBiol)$ProteinID == ProteinID)
	g <- as.numeric(resid_genes_TechBiol_sub3023[index.g,3:3025])
	p <- as.numeric(resid_proteins_TechBiol[index.p,3:3025])
	one_trio <- data.frame(g,p,age)
	colnames(one_trio)[1] <- 'Exprs'
	colnames(one_trio)[2] <- 'Protein'
	colnames(one_trio)[3] <- 'Age'
	s1 <- score(modelGENE, data = one_trio, type ='bic-g')
    s2 <- score(modelPROTEIN, data = one_trio, type ='bic-g')
    s3 <- score(modelINDEP, data = one_trio, type ='bic-g')
    BIC_Exprs_Prot[i,] <- c(GeneID,ProteinID,s1,s2,s3)
}

BIC_Exprs_Prot <- data.frame(BIC_Exprs_Prot)
for(i in 1:nrow(BIC_Exprs_Prot)){
	s1 <- as.numeric(BIC_Exprs_Prot[i,"modelGENE"])
	s2 <- as.numeric(BIC_Exprs_Prot[i,"modelPROTEIN"])
	s3 <- as.numeric(BIC_Exprs_Prot[i,"modelINDEP"])
	if(abs(s1-s2) > (s1-s3)){ # s2 < s3
		if(s3 < s1){	# s2 < s3 < s1
			BIC_Exprs_Prot$TOP[i] <- "modelPROTEIN" 
			BIC_Exprs_Prot$DIFF[i] <- abs(s2-s3)}	
		else if(s1 < s2){	# s1 < s2 < s3
			BIC_Exprs_Prot$TOP[i] <- "modelGENE" 
			BIC_Exprs_Prot$DIFF[i] <- abs(s1-s2)}
		else{ # s2 < s1 < s3
			BIC_Exprs_Prot$TOP[i] <- "modelPROTEIN" 
			BIC_Exprs_Prot$DIFF[i] <- abs(s2-s1)}}
	else{									# s3 < s2
		if(s2 < s1){	# s3 < s2 < s1
			BIC_Exprs_Prot$TOP[i] <- "modelINDEP" 
			BIC_Exprs_Prot$DIFF[i] <- abs(s3-s2)}
		else if(s1 < s3){ 	# s1 < s3 < s2
			BIC_Exprs_Prot$TOP[i] <- "modelGENE" 
			BIC_Exprs_Prot$DIFF[i] <- abs(s1-s3)}
		else{		# s3 < s1 < s2
			BIC_Exprs_Prot$TOP[i] <- "modelINDEP" 
			BIC_Exprs_Prot$DIFF[i] <- abs(s3-s1)}}
}
# The difference between the BICs for two models (∆BIC) can be computed,
# where guidelines indicate that a difference greater than 10 indicates very strong evidence;
# 6–10 indicates strong evidence;
# 2–6 indicates positive evidence;
# and 0–2 indicates weak evidence for the more complex model (Raftery, 1995)
for(i in 1:nrow(BIC_Exprs_Prot)){
	if(BIC_Exprs_Prot$DIFF[i] > 10){
		BIC_Exprs_Prot$SIGNIF[i]  <- "***" }
	else if(BIC_Exprs_Prot$DIFF[i]>6 && BIC_Exprs_Prot$DIFF[i]<10){
	BIC_Exprs_Prot$SIGNIF[i]  <- "**" }
	else if(BIC_Exprs_Prot$DIFF[i]>2 && BIC_Exprs_Prot$DIFF[i] <6){
	BIC_Exprs_Prot$SIGNIF[i]  <- "*" }
	else {BIC_Exprs_Prot$SIGNIF[i]  <- "-" }	
}
write.table(BIC_Exprs_Prot, file = paste("/users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Prot", index, sep="_"), quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
#colnames(BIC_Exprs_Prot) <- c("GeneID","ProteinID","modelGENE", "modelPROTEIN","modelINDEP","TOP","DIFF","SIGNIF")