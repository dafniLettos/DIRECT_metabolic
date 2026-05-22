#!/usr/bin/Rscript
index <- as.integer(Sys.getenv("PBS_ARRAYID"))
library("bnlearn")

#	Load Data
relat_TarMetTarMet <- read.table("/users/home/dafmic/CrossSectional//RELATIONS/results_TOTALS/TarMetTarMet_Assoc_Signif.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
relat_TarMetTarMet$MetaboliteID2 <- trimws(relat_TarMetTarMet$MetaboliteID2)
covar_Targf <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Metabolomics/Metabolites_Targeted_COVARIATES.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_tarMet_TechBiol <- read.table("/users/home/Teams/teamVIP/Analysis/Dafni/Baseline/Residuals/TarMet_TechBiol_withAgeEffects.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
resid_tarMet_TechBiol$Metabolite <- trimws(resid_tarMet_TechBiol$Metabolite)
#*******************************************************************************************		TARMET --- TARMET		******************************************************************************************************
modelINDEP <- empty.graph(c("Age", "Metab1", "Metab2"))
modelMETAB1 <- empty.graph(c("Age", "Metab1", "Metab2"))
modelMETAB2 <- empty.graph(c("Age", "Metab1", "Metab2"))
modelstring(modelINDEP) <- "[Metab1][Metab2][Age][Metab1|Age][Metab2|Age]" # Independet model
modelstring(modelMETAB1) <- "[Metab1][Metab2][Age][Metab1|Age][Metab2|Metab1]" # Metab1 via Metab2
modelstring(modelMETAB2) <- "[Metab1][Metab2][Age][Metab2|Age][Metab1|Metab2]" # Metab2 via Metab1
#	13,186
BIC_TarMet_TarMet <- matrix(0,ncol=5,nrow=1319)
colnames(BIC_TarMet_TarMet) <- c("MetaboliteID1","MetaboliteID2","modelMETAB1", "modelMETAB2","modelINDEP")
relat_TarMetTarMet2 <- relat_TarMetTarMet[(index*1319-1318):min(index*1319,nrow(relat_TarMetTarMet)),]

for(i in 1:nrow(relat_TarMetTarMet2)){
	MetabID1 <- as.character(relat_TarMetTarMet2[i,1])
	MetabID2 <- as.character(relat_TarMetTarMet2[i,2])
	age <- as.numeric(covar_Targf$Age)
	index.g1 <- which(resid_tarMet_TechBiol$Metabolite == MetabID1)
	index.g2 <- which(resid_tarMet_TechBiol$Metabolite == MetabID2)
	g1 <- as.numeric(resid_tarMet_TechBiol[index.g1,2:3028])
	g2 <- as.numeric(resid_tarMet_TechBiol[index.g2,2:3028])
	one_trio <- data.frame(g1,g2,age)
	colnames(one_trio)[1] <- 'Metab1'
	colnames(one_trio)[2] <- 'Metab2'
	colnames(one_trio)[3] <- 'Age'
	s1 <- score(modelMETAB1, data = one_trio, type ='bic-g')
    s2 <- score(modelMETAB2, data = one_trio, type ='bic-g')
    s3 <- score(modelINDEP, data = one_trio, type ='bic-g')
    BIC_TarMet_TarMet[i,] <- c(MetabID1,MetabID2,s1,s2,s3)
}

BIC_TarMet_TarMet <- data.frame(BIC_TarMet_TarMet)
for(i in 1:nrow(BIC_TarMet_TarMet)){
	s1 <- as.numeric(BIC_TarMet_TarMet[i,"modelMETAB1"])
	s2 <- as.numeric(BIC_TarMet_TarMet[i,"modelMETAB2"])
	s3 <- as.numeric(BIC_TarMet_TarMet[i,"modelINDEP"])
	if(BIC_TarMet_TarMet[i,1] == BIC_TarMet_TarMet[i,2]){
		BIC_TarMet_TarMet$TOP[i] <- "1-1"
		BIC_TarMet_TarMet$DIFF[i] <- 1}
	else if(abs(s1-s2) > (s1-s3)){ # s2 < s3
		if(s3 < s1){	# s2 < s3 < s1
			BIC_TarMet_TarMet$TOP[i] <- "modelMETAB2" 
			BIC_TarMet_TarMet$DIFF[i] <- abs(s2-s3)}	
		else if(s1 < s2){	# s1 < s2 < s3
			BIC_TarMet_TarMet$TOP[i] <- "modelMETAB1" 
			BIC_TarMet_TarMet$DIFF[i] <- abs(s1-s2)}
		else{ # s2 < s1 < s3
			BIC_TarMet_TarMet$TOP[i] <- "modelMETAB2" 
			BIC_TarMet_TarMet$DIFF[i] <- abs(s2-s1)}}
	else{									# s3 < s2
		if(s2 < s1){	# s3 < s2 < s1
			BIC_TarMet_TarMet$TOP[i] <- "modelINDEP" 
			BIC_TarMet_TarMet$DIFF[i] <- abs(s3-s2)}
		else if(s1 < s3){ 	# s1 < s3 < s2
			BIC_TarMet_TarMet$TOP[i] <- "modelMETAB1" 
			BIC_TarMet_TarMet$DIFF[i] <- abs(s1-s3)}
		else{		# s3 < s1 < s2
			BIC_TarMet_TarMet$TOP[i] <- "modelINDEP" 
			BIC_TarMet_TarMet$DIFF[i] <- abs(s3-s1)}}
}
# The difference between the BICs for two models (∆BIC) can be computed,
# where guidelines indicate that a difference greater than 10 indicates very strong evidence;
# 6–10 indicates strong evidence;
# 2–6 indicates positive evidence;
# and 0–2 indicates weak evidence for the more complex model (Raftery, 1995)
for(i in 1:nrow(BIC_TarMet_TarMet)){
	if(BIC_TarMet_TarMet$DIFF[i] > 10){
		BIC_TarMet_TarMet$SIGNIF[i]  <- "***" }
	else if(BIC_TarMet_TarMet$DIFF[i]>6 && BIC_TarMet_TarMet$DIFF[i]<10){
	BIC_TarMet_TarMet$SIGNIF[i]  <- "**" }
	else if(BIC_TarMet_TarMet$DIFF[i]>2 && BIC_TarMet_TarMet$DIFF[i] <6){
	BIC_TarMet_TarMet$SIGNIF[i]  <- "*" }
	else {BIC_TarMet_TarMet$SIGNIF[i]  <- "-" }	
}
write.table(BIC_TarMet_TarMet, file = paste("/users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_TarMet_TarMet", index, sep="_"), quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)
#colnames(BIC_TarMet_TarMet) <- c("MetaboliteID1","MetaboliteID2","modelMETAB1", "modelMETAB2","modelINDEP","TOP","DIFF","SIGNIF")