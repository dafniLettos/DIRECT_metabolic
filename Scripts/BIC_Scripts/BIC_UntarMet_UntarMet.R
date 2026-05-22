#!/usr/bin/Rscript
index <- as.integer(Sys.getenv("PBS_ARRAYID"))
library("bnlearn")

#	Load Data
relat_UntarMetUntarMet <- read.table("/users/home/dafmic/CrossSectional//RELATIONS/results_TOTALS/UntarMetUntarMet_Assoc_Signif.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
relat_UntarMetUntarMet$MetaboliteID2 <- trimws(relat_UntarMetUntarMet$MetaboliteID2)
covar_Untargf <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Metabolomics/Metabolites_Untargeted_COVARIATES.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_untarMet_TechBiol <- read.delim("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/UntarMet_Residuals_TechBiol_withAgeEffects.txt",header = TRUE)
resid_untarMet_TechBiol$BIOCHEMICAL2 <- trimws(resid_untarMet_TechBiol$BIOCHEMICAL2)
#*******************************************************************************************		TARMET --- TARMET		******************************************************************************************************
modelINDEP <- empty.graph(c("Age", "Metab1", "Metab2"))
modelMETAB1 <- empty.graph(c("Age", "Metab1", "Metab2"))
modelMETAB2 <- empty.graph(c("Age", "Metab1", "Metab2"))
modelstring(modelINDEP) <- "[Metab1][Metab2][Age][Metab1|Age][Metab2|Age]" # Independet model
modelstring(modelMETAB1) <- "[Metab1][Metab2][Age][Metab1|Age][Metab2|Metab1]" # Metab1 via Metab2
modelstring(modelMETAB2) <- "[Metab1][Metab2][Age][Metab2|Age][Metab1|Metab2]" # Metab2 via Metab1
#	7,617
BIC_UntarMet_UntarMet <- matrix(0,ncol=5,nrow=762)
colnames(BIC_UntarMet_UntarMet) <- c("MetaboliteID1","MetaboliteID2","modelMETAB1", "modelMETAB2","modelINDEP")
relat_UntarMetUntarMet2 <- relat_UntarMetUntarMet[(index*762-761):min(index*762,nrow(relat_UntarMetUntarMet)),]

for(i in 1:nrow(relat_UntarMetUntarMet2)){
	MetabID1 <- as.character(relat_UntarMetUntarMet2[i,1])
	MetabID2 <- as.character(relat_UntarMetUntarMet2[i,3])
	age <- as.numeric(covar_Untargf$Age)
	index.g1 <- which(resid_untarMet_TechBiol$BIOCHEMICAL2 == MetabID1)
	index.g2 <- which(resid_untarMet_TechBiol$BIOCHEMICAL2 == MetabID2)
	g1 <- as.numeric(resid_untarMet_TechBiol[index.g1,3:2998])
	g2 <- as.numeric(resid_untarMet_TechBiol[index.g2,3:2998])
	one_trio <- data.frame(g1,g2,age)
	colnames(one_trio)[1] <- 'Metab1'
	colnames(one_trio)[2] <- 'Metab2'
	colnames(one_trio)[3] <- 'Age'
	s1 <- score(modelMETAB1, data = one_trio, type ='bic-g')
    s2 <- score(modelMETAB2, data = one_trio, type ='bic-g')
    s3 <- score(modelINDEP, data = one_trio, type ='bic-g')
    BIC_UntarMet_UntarMet[i,] <- c(MetabID1,MetabID2,s1,s2,s3)
}

BIC_UntarMet_UntarMet <- data.frame(BIC_UntarMet_UntarMet)
for(i in 1:nrow(BIC_UntarMet_UntarMet)){
	s1 <- as.numeric(BIC_UntarMet_UntarMet[i,"modelMETAB1"])
	s2 <- as.numeric(BIC_UntarMet_UntarMet[i,"modelMETAB2"])
	s3 <- as.numeric(BIC_UntarMet_UntarMet[i,"modelINDEP"])
	if(BIC_UntarMet_UntarMet[i,1] == BIC_UntarMet_UntarMet[i,2]){
		BIC_UntarMet_UntarMet$TOP[i] <- "1-1"
		BIC_UntarMet_UntarMet$DIFF[i] <- 1}
	else if(abs(s1-s2) > (s1-s3)){ # s2 < s3
		if(s3 < s1){	# s2 < s3 < s1
			BIC_UntarMet_UntarMet$TOP[i] <- "modelMETAB2" 
			BIC_UntarMet_UntarMet$DIFF[i] <- abs(s2-s3)}	
		else if(s1 < s2){	# s1 < s2 < s3
			BIC_UntarMet_UntarMet$TOP[i] <- "modelMETAB1" 
			BIC_UntarMet_UntarMet$DIFF[i] <- abs(s1-s2)}
		else{ # s2 < s1 < s3
			BIC_UntarMet_UntarMet$TOP[i] <- "modelMETAB2" 
			BIC_UntarMet_UntarMet$DIFF[i] <- abs(s2-s1)}}
	else{									# s3 < s2
		if(s2 < s1){	# s3 < s2 < s1
			BIC_UntarMet_UntarMet$TOP[i] <- "modelINDEP" 
			BIC_UntarMet_UntarMet$DIFF[i] <- abs(s3-s2)}
		else if(s1 < s3){ 	# s1 < s3 < s2
			BIC_UntarMet_UntarMet$TOP[i] <- "modelMETAB1" 
			BIC_UntarMet_UntarMet$DIFF[i] <- abs(s1-s3)}
		else{		# s3 < s1 < s2
			BIC_UntarMet_UntarMet$TOP[i] <- "modelINDEP" 
			BIC_UntarMet_UntarMet$DIFF[i] <- abs(s3-s1)}}
}
# The difference between the BICs for two models (∆BIC) can be computed,
# where guidelines indicate that a difference greater than 10 indicates very strong evidence;
# 6–10 indicates strong evidence;
# 2–6 indicates positive evidence;
# and 0–2 indicates weak evidence for the more complex model (Raftery, 1995)
for(i in 1:nrow(BIC_UntarMet_UntarMet)){
	if(BIC_UntarMet_UntarMet$DIFF[i] > 10){
		BIC_UntarMet_UntarMet$SIGNIF[i]  <- "***" }
	else if(BIC_UntarMet_UntarMet$DIFF[i]>6 && BIC_UntarMet_UntarMet$DIFF[i]<10){
	BIC_UntarMet_UntarMet$SIGNIF[i]  <- "**" }
	else if(BIC_UntarMet_UntarMet$DIFF[i]>2 && BIC_UntarMet_UntarMet$DIFF[i] <6){
	BIC_UntarMet_UntarMet$SIGNIF[i]  <- "*" }
	else {BIC_UntarMet_UntarMet$SIGNIF[i]  <- "-" }	
}
write.table(BIC_UntarMet_UntarMet, file = paste("/users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_UntarMet_UntarMet", index, sep="_"), quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
#colnames(BIC_UntarMet_UntarMet) <- c("MetaboliteID1","MetaboliteID2","modelMETAB1", "modelMETAB2","modelINDEP","TOP","DIFF","SIGNIF")