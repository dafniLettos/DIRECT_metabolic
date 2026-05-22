#!/usr/bin/Rscript
index <- as.integer(Sys.getenv("PBS_ARRAYID"))
library("bnlearn")

#	Load Data
relat_ExprsTarMet <- read.table("/users/home/dafmic/CrossSectional//RELATIONS/results_TOTALS/ExprsTarMet_Assoc_Signif.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
relat_ExprsTarMet$GeneID <- trimws(relat_ExprsTarMet$GeneID)
relat_ExprsTarMet$MetaboliteID <- trimws(relat_ExprsTarMet$MetaboliteID)
covar_Exprsf <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Transcriptomics/GeneExpression_COVARIATES.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_genes_TechBiol <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/Exprs_Residuals_TechBiol_withAgeEffects.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_genes_TechBiol$GeneID <- trimws(resid_genes_TechBiol$GeneID)
resid_tarMet_TechBiol <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/TarMet_TechBiol_withAgeEffects.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_tarMet_TechBiol$Metabolite <- trimws(resid_tarMet_TechBiol$Metabolite)
#*******************************************************************************************		EXPRS  -- TARMET		******************************************************************************************************
modelINDEP <- empty.graph(c("Age", "Exprs", "TarMet"))
modelTARMET <- empty.graph(c("Age", "Exprs", "TarMet"))
modelGENE <- empty.graph(c("Age", "Exprs", "TarMet"))
modelstring(modelINDEP) <- "[Exprs][TarMet][Age][Exprs|Age][TarMet|Age]" # Independet model
modelstring(modelTARMET) <- "[Exprs][TarMet][Age][Exprs|Age][TarMet|Exprs]" # Gene via TarMet
modelstring(modelGENE) <- "[Exprs][TarMet][Age][TarMet|Age][Exprs|TarMet]" # TarMet via Gene
#	217,970
BIC_Exprs_TarMet <- matrix(0,ncol=5,nrow=2180)
colnames(BIC_Exprs_TarMet) <- c("GeneID","MetaboliteID","modelGENE", "modelTARMET","modelINDEP")
relat_ExprsTarMet2 <- relat_ExprsTarMet[(index*2180-2179):min(index*2180,nrow(relat_ExprsTarMet)),]

for(i in 1:nrow(relat_ExprsTarMet2)){
	GeneID <- as.character(relat_ExprsTarMet2[i,1])
	MetaboliteID <- as.character(relat_ExprsTarMet2[i,3])
	age <- as.numeric(covar_Exprsf$Age)
	index.g <- which(as.data.frame(resid_genes_TechBiol)$GeneID == GeneID)
	index.m <- which(as.data.frame(resid_tarMet_TechBiol)$Metabolite == MetaboliteID)
	g <- as.numeric(resid_genes_TechBiol[index.g,3:3029])
	m <- as.numeric(resid_tarMet_TechBiol[index.m,2:3028])
	one_trio <- data.frame(g,m,age)
	colnames(one_trio)[1] <- 'Exprs'
	colnames(one_trio)[2] <- 'TarMet'
	colnames(one_trio)[3] <- 'Age'
	s1 <- score(modelGENE, data = one_trio, type ='bic-g')
    s2 <- score(modelTARMET, data = one_trio, type ='bic-g')
    s3 <- score(modelINDEP, data = one_trio, type ='bic-g')
    BIC_Exprs_TarMet[i,] <- c(GeneID,MetaboliteID,s1,s2,s3)
}

BIC_Exprs_TarMet <- data.frame(BIC_Exprs_TarMet)
for(i in 1:nrow(BIC_Exprs_TarMet)){
	s1 <- as.numeric(BIC_Exprs_TarMet[i,"modelGENE"])
	s2 <- as.numeric(BIC_Exprs_TarMet[i,"modelTARMET"])
	s3 <- as.numeric(BIC_Exprs_TarMet[i,"modelINDEP"])
	if(abs(s1-s2) > (s1-s3)){ # s2 < s3
		if(s3 < s1){	# s2 < s3 < s1
			BIC_Exprs_TarMet$TOP[i] <- "modelTARMET" 
			BIC_Exprs_TarMet$DIFF[i] <- abs(s2-s3)}	
		else if(s1 < s2){	# s1 < s2 < s3
			BIC_Exprs_TarMet$TOP[i] <- "modelGENE" 
			BIC_Exprs_TarMet$DIFF[i] <- abs(s1-s2)}
		else{ # s2 < s1 < s3
			BIC_Exprs_TarMet$TOP[i] <- "modelTARMET" 
			BIC_Exprs_TarMet$DIFF[i] <- abs(s2-s1)}}
	else{									# s3 < s2
		if(s2 < s1){	# s3 < s2 < s1
			BIC_Exprs_TarMet$TOP[i] <- "modelINDEP" 
			BIC_Exprs_TarMet$DIFF[i] <- abs(s3-s2)}
		else if(s1 < s3){ 	# s1 < s3 < s2
			BIC_Exprs_TarMet$TOP[i] <- "modelGENE" 
			BIC_Exprs_TarMet$DIFF[i] <- abs(s1-s3)}
		else{		# s3 < s1 < s2
			BIC_Exprs_TarMet$TOP[i] <- "modelINDEP" 
			BIC_Exprs_TarMet$DIFF[i] <- abs(s3-s1)}}
}
# The difference between the BICs for two models (∆BIC) can be computed,
# where guidelines indicate that a difference greater than 10 indicates very strong evidence;
# 6–10 indicates strong evidence;
# 2–6 indicates positive evidence;
# and 0–2 indicates weak evidence for the more complex model (Raftery, 1995)
for(i in 1:nrow(BIC_Exprs_TarMet)){
	if(BIC_Exprs_TarMet$DIFF[i] > 10){
		BIC_Exprs_TarMet$SIGNIF[i]  <- "***" }
	else if(BIC_Exprs_TarMet$DIFF[i]>6 && BIC_Exprs_TarMet$DIFF[i]<10){
	BIC_Exprs_TarMet$SIGNIF[i]  <- "**" }
	else if(BIC_Exprs_TarMet$DIFF[i]>2 && BIC_Exprs_TarMet$DIFF[i] <6){
	BIC_Exprs_TarMet$SIGNIF[i]  <- "*" }
	else {BIC_Exprs_TarMet$SIGNIF[i]  <- "-" }	
}

write.table(BIC_Exprs_TarMet, file = paste("/users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_TarMet", index, sep="_"), quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
#colnames(BIC_Exprs_TarMet) <- c("GeneID","MetaboliteID","modelGENE", "modelTARMET","modelINDEP",TOP","DIFF","SIGNIF")