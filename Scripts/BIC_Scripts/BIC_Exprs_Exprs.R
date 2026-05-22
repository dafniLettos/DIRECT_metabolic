#!/usr/bin/Rscript
index <- as.integer(Sys.getenv("PBS_ARRAYID"))
library("bnlearn")

#	Load Data
relat_ExprsExprs <- read.table("/users/home/dafmic/CrossSectional//RELATIONS/results_TOTALS/ExprsExprs_Assoc_Signif_aa",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
#relat_ExprsExprs$GeneID1 <- trimws(relat_ExprsExprs$GeneID1)
relat_ExprsExprs$GeneID2 <- trimws(relat_ExprsExprs$GeneID2)
covar_Exprsf <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Transcriptomics/GeneExpression_COVARIATES.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_genes_TechBiol <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/Exprs_Residuals_TechBiol_withAgeEffects.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_genes_TechBiol$GeneID <- trimws(resid_genes_TechBiol$GeneID)
#*******************************************************************************************		EXPRS --- EXPRS		******************************************************************************************************
modelINDEP <- empty.graph(c("Age", "Exprs1", "Exprs2"))
modelGENE1 <- empty.graph(c("Age", "Exprs1", "Exprs2"))
modelGENE2 <- empty.graph(c("Age", "Exprs1", "Exprs2"))
modelstring(modelINDEP) <- "[Exprs1][Exprs2][Age][Exprs1|Age][Exprs2|Age]" # Independet model
modelstring(modelGENE1) <- "[Exprs1][Exprs2][Age][Exprs1|Age][Exprs2|Exprs1]" # Gene1 via Gene2
modelstring(modelGENE2) <- "[Exprs1][Exprs2][Age][Exprs2|Age][Exprs1|Exprs2]" # Gene2 via Gene1
#	1,542,867 / 108,000,652
BIC_Exprs_Exprs <- matrix(0,ncol=5,nrow=1543)
colnames(BIC_Exprs_Exprs) <- c("GeneID1","GeneID2","modelGENE1", "modelGENE2","modelINDEP")
relat_ExprsExprs2 <- relat_ExprsExprs[(index*1543-1542):min(index*1543,nrow(relat_ExprsExprs)),]

for(i in 1:nrow(relat_ExprsExprs2)){
	GeneID1 <- as.character(relat_ExprsExprs2[i,1])
	GeneID2 <- as.character(relat_ExprsExprs2[i,3])
	age <- as.numeric(covar_Exprsf$Age)
	index.g1 <- which(resid_genes_TechBiol$GeneID == GeneID1)
	index.g2 <- which(resid_genes_TechBiol$GeneID == GeneID2)
	g1 <- as.numeric(resid_genes_TechBiol[index.g1,3:3029])
	g2 <- as.numeric(resid_genes_TechBiol[index.g2,3:3029])
	one_trio <- data.frame(g1,g2,age)
	colnames(one_trio)[1] <- 'Exprs1'
	colnames(one_trio)[2] <- 'Exprs2'
	colnames(one_trio)[3] <- 'Age'
	s1 <- score(modelGENE1, data = one_trio, type ='bic-g')
    s2 <- score(modelGENE2, data = one_trio, type ='bic-g')
    s3 <- score(modelINDEP, data = one_trio, type ='bic-g')
    BIC_Exprs_Exprs[i,] <- c(GeneID1,GeneID2,s1,s2,s3)
}

BIC_Exprs_Exprs <- data.frame(BIC_Exprs_Exprs)
for(i in 1:nrow(BIC_Exprs_Exprs)){
	s1 <- as.numeric(BIC_Exprs_Exprs[i,"modelGENE1"])
	s2 <- as.numeric(BIC_Exprs_Exprs[i,"modelGENE2"])
	s3 <- as.numeric(BIC_Exprs_Exprs[i,"modelINDEP"])
	if(BIC_Exprs_Exprs[i,1] == BIC_Exprs_Exprs[i,2]){
		BIC_Exprs_Exprs$TOP[i] <- "1-1"
		BIC_Exprs_Exprs$DIFF[i] <- 1}
	else if(abs(s1-s2) > (s1-s3)){ # s2 < s3
		if(s3 < s1){	# s2 < s3 < s1
			BIC_Exprs_Exprs$TOP[i] <- "modelGENE2" 
			BIC_Exprs_Exprs$DIFF[i] <- abs(s2-s3)}	
		else if(s1 < s2){	# s1 < s2 < s3
			BIC_Exprs_Exprs$TOP[i] <- "modelGENE1" 
			BIC_Exprs_Exprs$DIFF[i] <- abs(s1-s2)}
		else{ # s2 < s1 < s3
			BIC_Exprs_Exprs$TOP[i] <- "modelGENE2" 
			BIC_Exprs_Exprs$DIFF[i] <- abs(s2-s1)}}
	else{									# s3 < s2
		if(s2 < s1){	# s3 < s2 < s1
			BIC_Exprs_Exprs$TOP[i] <- "modelINDEP" 
			BIC_Exprs_Exprs$DIFF[i] <- abs(s3-s2)}
		else if(s1 < s3){ 	# s1 < s3 < s2
			BIC_Exprs_Exprs$TOP[i] <- "modelGENE1" 
			BIC_Exprs_Exprs$DIFF[i] <- abs(s1-s3)}
		else{		# s3 < s1 < s2
			BIC_Exprs_Exprs$TOP[i] <- "modelINDEP" 
			BIC_Exprs_Exprs$DIFF[i] <- abs(s3-s1)}}
}
# The difference between the BICs for two models (∆BIC) can be computed,
# where guidelines indicate that a difference greater than 10 indicates very strong evidence;
# 6–10 indicates strong evidence;
# 2–6 indicates positive evidence;
# and 0–2 indicates weak evidence for the more complex model (Raftery, 1995)
for(i in 1:nrow(BIC_Exprs_Exprs)){
	if(BIC_Exprs_Exprs$DIFF[i] > 10){
		BIC_Exprs_Exprs$SIGNIF[i]  <- "***" }
	else if(BIC_Exprs_Exprs$DIFF[i]>6 && BIC_Exprs_Exprs$DIFF[i]<10){
	BIC_Exprs_Exprs$SIGNIF[i]  <- "**" }
	else if(BIC_Exprs_Exprs$DIFF[i]>2 && BIC_Exprs_Exprs$DIFF[i] <6){
	BIC_Exprs_Exprs$SIGNIF[i]  <- "*" }
	else {BIC_Exprs_Exprs$SIGNIF[i]  <- "-" }	
}
write.table(BIC_Exprs_Exprs, file = paste("/users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_aa", index, sep="_"), quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
#colnames(BIC_Exprs_Exprs) <- c("GeneID1","GeneID2","modelGENE1", "modelGENE2","modelINDEP","TOP","DIFF","SIGNIF")