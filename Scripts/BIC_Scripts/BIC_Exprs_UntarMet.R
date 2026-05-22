#!/usr/bin/Rscript
index <- as.integer(Sys.getenv("PBS_ARRAYID"))
library("bnlearn")

#	Load Data
raw_genesE <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Transcriptomics/GeneExpression_processed.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
relat_ExprsUntarMet <- read.table("/users/home/dafmic/CrossSectional//RELATIONS/results_TOTALS/ExprsUntarMet_Assoc_Signif.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
relat_ExprsUntarMet$GeneID <- trimws(relat_ExprsUntarMet$GeneID )
relat_ExprsUntarMet$BIOCHEMICAL2 <- trimws(relat_ExprsUntarMet$BIOCHEMICAL2)
raw_met_untarg_Imp <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Metabolomics/Metabolites_Untargeted_processed.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
covar_Untargf <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Metabolomics/Metabolites_Untargeted_COVARIATES.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
rownames(covar_Untargf) <- colnames(raw_met_untarg_Imp)
v <- rownames(covar_Untargf)[order(match(rownames(covar_Untargf), colnames(raw_genesE)))]
covar_Untargf <- covar_Untargf[match(v,rownames(covar_Untargf)),]
raw_met_untarg_Imp <- raw_met_untarg_Imp[,rownames(covar_Untargf)]
resid_genes_TechBiol <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/Exprs_Residuals_TechBiol_withAgeEffects.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_genes_TechBiol$GeneID <- trimws(resid_genes_TechBiol$GeneID)
resid_untarMet_TechBiol <- read.delim("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/UntarMet_Residuals_TechBiol_withAgeEffects.txt", header = TRUE)
resid_untarMet_TechBiol$BIOCHEMICAL2 <- trimws(resid_untarMet_TechBiol$BIOCHEMICAL2)
resid_genes_TechBiol_sub2996 <- data.frame(resid_genes_TechBiol[,1:2],resid_genes_TechBiol[,c(which(colnames(resid_genes_TechBiol) %in% colnames(raw_met_untarg_Imp)))])
colnames(resid_genes_TechBiol_sub2996)[3:2998] <- colnames(raw_met_untarg_Imp)
#*******************************************************************************************		EXPRS --- UNTARMET	******************************************************************************************************
modelINDEP <- empty.graph(c("Age", "Exprs", "UntarMet"))
modelUNTARMET <- empty.graph(c("Age", "Exprs", "UntarMet"))
modelGENE <- empty.graph(c("Age", "Exprs", "UntarMet"))
modelstring(modelINDEP) <- "[Exprs][UntarMet][Age][Exprs|Age][UntarMet|Age]" # Independet model
modelstring(modelUNTARMET) <- "[Exprs][UntarMet][Age][Exprs|Age][UntarMet|Exprs]" # Gene via UntarMet
modelstring(modelGENE) <- "[Exprs][UntarMet][Age][UntarMet|Age][Exprs|UntarMet]" # UntarMet via Gene
#	111,944
BIC_Exprs_UntarMet <- matrix(0,ncol=5,nrow=1120)
colnames(BIC_Exprs_UntarMet) <- c("GeneID","MetaboliteID","modelGENE", "modelUNTARMET","modelINDEP")
relat_ExprsUntarMet2 <- relat_ExprsUntarMet[(index*1120-1119):min(index*1120,nrow(relat_ExprsUntarMet)),]

for(i in 1:nrow(relat_ExprsUntarMet2)){
	GeneID <- as.character(relat_ExprsUntarMet2[i,1])
	MetaboliteID <- as.character(relat_ExprsUntarMet2[i,3])
	age <- as.numeric(covar_Untargf$Age)
	index.g <- which(as.data.frame(resid_genes_TechBiol_sub2996)$GeneID == GeneID)
	index.m <- which(as.data.frame(resid_untarMet_TechBiol)$BIOCHEMICAL2 == MetaboliteID)
	g <- as.numeric(resid_genes_TechBiol_sub2996[index.g,3:2998])
	m <- as.numeric(resid_untarMet_TechBiol[index.m,3:2998])
	one_trio <- data.frame(g,m,age)
	colnames(one_trio)[1] <- 'Exprs'
	colnames(one_trio)[2] <- 'UntarMet'
	colnames(one_trio)[3] <- 'Age'
	s1 <- score(modelGENE, data = one_trio, type ='bic-g')
    s2 <- score(modelUNTARMET, data = one_trio, type ='bic-g')
    s3 <- score(modelINDEP, data = one_trio, type ='bic-g')
    BIC_Exprs_UntarMet[i,] <- c(GeneID,MetaboliteID,s1,s2,s3)
}

BIC_Exprs_UntarMet <- data.frame(BIC_Exprs_UntarMet)
for(i in 1:nrow(BIC_Exprs_UntarMet)){
	s1 <- as.numeric(BIC_Exprs_UntarMet[i,"modelGENE"])
	s2 <- as.numeric(BIC_Exprs_UntarMet[i,"modelUNTARMET"])
	s3 <- as.numeric(BIC_Exprs_UntarMet[i,"modelINDEP"])
	if(abs(s1-s2) > (s1-s3)){ # s2 < s3
		if(s3 < s1){	# s2 < s3 < s1
			BIC_Exprs_UntarMet$TOP[i] <- "modelUNTARMET" 
			BIC_Exprs_UntarMet$DIFF[i] <- abs(s2-s3)}	
		else if(s1 < s2){	# s1 < s2 < s3
			BIC_Exprs_UntarMet$TOP[i] <- "modelGENE" 
			BIC_Exprs_UntarMet$DIFF[i] <- abs(s1-s2)}
		else{ # s2 < s1 < s3
			BIC_Exprs_UntarMet$TOP[i] <- "modelUNTARMET" 
			BIC_Exprs_UntarMet$DIFF[i] <- abs(s2-s1)}}
	else{									# s3 < s2
		if(s2 < s1){	# s3 < s2 < s1
			BIC_Exprs_UntarMet$TOP[i] <- "modelINDEP" 
			BIC_Exprs_UntarMet$DIFF[i] <- abs(s3-s2)}
		else if(s1 < s3){ 	# s1 < s3 < s2
			BIC_Exprs_UntarMet$TOP[i] <- "modelGENE" 
			BIC_Exprs_UntarMet$DIFF[i] <- abs(s1-s3)}
		else{		# s3 < s1 < s2
			BIC_Exprs_UntarMet$TOP[i] <- "modelINDEP" 
			BIC_Exprs_UntarMet$DIFF[i] <- abs(s3-s1)}}
}
# The difference between the BICs for two models (∆BIC) can be computed,
# where guidelines indicate that a difference greater than 10 indicates very strong evidence;
# 6–10 indicates strong evidence;
# 2–6 indicates positive evidence;
# and 0–2 indicates weak evidence for the more complex model (Raftery, 1995)
for(i in 1:nrow(BIC_Exprs_UntarMet)){
	if(BIC_Exprs_UntarMet$DIFF[i] > 10){
		BIC_Exprs_UntarMet$SIGNIF[i]  <- "***" }
	else if(BIC_Exprs_UntarMet$DIFF[i]>6 && BIC_Exprs_UntarMet$DIFF[i]<10){
	BIC_Exprs_UntarMet$SIGNIF[i]  <- "**" }
	else if(BIC_Exprs_UntarMet$DIFF[i]>2 && BIC_Exprs_UntarMet$DIFF[i] <6){
	BIC_Exprs_UntarMet$SIGNIF[i]  <- "*" }
	else {BIC_Exprs_UntarMet$SIGNIF[i]  <- "-" }	
}

write.table(BIC_Exprs_UntarMet, file = paste("/users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_UntarMet", index, sep="_"), quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
#colnames(BIC_Exprs_UntarMet) <- c("GeneID","MetaboliteID","modelGENE", "modelUNTARMET","modelINDEP","TOP","DIFF","SIGNIF")