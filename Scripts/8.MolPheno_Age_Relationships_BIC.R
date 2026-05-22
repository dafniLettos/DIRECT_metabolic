# THIS FILE CONTAINS THE R CODE USED TO IDENTIFY RELATIONSHIPS (ASSOCIATIONS) BETWEEN MOLECULAR PHENOTYPES 
# AND CONTEXT-DEPENDENT EFFECTS (OF AGE, SEX, BMI, HbA1c) MODULATING THESE RELATIONSHIPS
# FOR GENES, PROTEINS AND METABOLITES (TARGETED & UNTARGETED DATASET)
# USING LINEAR REGRESSION AND INTERACTION MODELS
# AND APPLYING MULTIPLE CORRECTION METHODS
# ALSO INCLUDES CODE FOR CONSTRUCTING BARPLOT AND UPSET PLOTS

#-------------------------------------------------------------------	LIBRARIES/PACKAGES USED	-------------------------------------------------------------------------
library(tibble)
library(dplyr)
library(qvalue)
library(data.table)
setDTthreads(threads = 0)
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Assoc_ExprsExprs <- fread("ExprsExpr_Assoc_AGE.txt")
Assoc_ExprsProt <- fread("ExprsProt_Assoc_AGE.txt")
Assoc_ExprsTarMet <- fread("ExprsTarMet_Assoc_AGE.txt")
Assoc_ProtProt <- fread("ProtProt_Assoc_AGE.txt")
Assoc_ExprsUntarMet <- fread("ExprsUntarMet_Assoc_AGE.txt")
Assoc_ProtTarMet <- fread("ProtTarMet_Assoc_AGE.txt")
Assoc_ProtUntarMet <- fread("ProtUntarMet_Assoc_AGE.txt")
Assoc_TarMetTarMet <- fread("TarMetTarMet_Assoc_AGE.txt")
Assoc_UntarMetUntarMet <- fread("UntarMetUntarMet_Assoc_AGE.txt")
Assoc_TarMetUntarMet <- fread("TarMetUntarMet_Assoc_AGE.txt")

Assoc_ExprsExprs$GeneNameA <- paste0(Assoc_ExprsExprs$GeneName1, "e")
Assoc_ExprsExprs$GeneNameB <- paste0(Assoc_ExprsExprs$GeneName2, "e")
Assoc_ExprsProt$GeneNameA <- paste0(Assoc_ExprsProt$GeneName, "e")
Assoc_ExprsProt$GeneNameB <- paste0(Assoc_ExprsProt$ProteinName, "p")
Assoc_ExprsTarMet$GeneNameA <- paste0(Assoc_ExprsTarMet$GeneName, "e")
Assoc_ExprsUntarMet$GeneNameA <- paste0(Assoc_ExprsUntarMet$GeneName, "e")
Assoc_ProtProt$GeneNameA <- paste0(Assoc_ProtProt$GeneName1, "p")
Assoc_ProtProt$GeneNameB <- paste0(Assoc_ProtProt$GeneName2, "p")
Assoc_ProtTarMet$GeneName2 <- paste0(Assoc_ProtTarMet$GeneName, "p")

Assoc_ExprsExprs <- Assoc_ExprsExprs[,c(9:10,5:8)]
colnames(Assoc_ExprsExprs) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_ExprsProt <- Assoc_ExprsProt[,c(10:11,5:8)]
colnames(TEMP_Assoc_ExprsProt) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_ProtProt <- Assoc_ProtProt[,c(10:11,5:8)]
colnames(TEMP_Assoc_ProtProt) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_ExprsTarMet <- Assoc_ExprsTarMet[,c(9,3:7)]
colnames(TEMP_Assoc_ExprsTarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_ExprsUntarMet <- Assoc_ExprsUntarMet[,c(10,3,5:8)]
colnames(TEMP_Assoc_ExprsUntarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_ProtTarMet <- Assoc_ProtTarMet[,c(9,3:7)]
colnames(TEMP_Assoc_ProtTarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_TarMetTarMet <- Assoc_TarMetTarMet[,c(1:6)]
colnames(TEMP_Assoc_TarMetTarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_UntarMetUntarMet <- Assoc_UntarMetUntarMet[,c(1,3,5:8)]
colnames(TEMP_Assoc_UntarMetUntarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_TarMetUntarMet <- Assoc_TarMetUntarMet[,c(1,2,4:7)]
colnames(TEMP_Assoc_TarMetUntarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
ALL <- rbind(Assoc_ExprsExprs, TEMP_Assoc_ExprsTarMet, TEMP_Assoc_ExprsUntarMet, TEMP_Assoc_ProtProt, TEMP_Assoc_ProtTarMet, TEMP_Assoc_TarMetTarMet, TEMP_Assoc_UntarMetUntarMet, TEMP_Assoc_TarMetUntarMet)

write.table(ALL, file = "S77_Relationships_AgeRelatedMolecules.tab", quote = FALSE, sep = "\t", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
Associations_ALL <- fread("S77_Relationships_AgeRelatedMolecules.tab")

#==================================================================
#                                           NETWORK
#==================================================================
setwd("~/Desktop/PhD/0.DIRECT_SEX_PAPER/NETWORK")
# GET NUMBERS OF SIGNIFICANT ASSOCIATIONS
sum(Assoc_ExprsExprs$Exprs_FDR < 0.05)
#	96,737,904/113,262,806  [85.4%]
sum(Assoc_ExprsProt$FDR < 0.05)
#	669,994/2,841,681 [23.6%]
sum(Assoc_ExprsTarMet$FDR < 0.05)
#	152,101/1,085,586 [14%]
sum(Assoc_ExprsUntarMet$FDR < 0.05)
#	151,747/2,000,884 [7.6%]
sum(Assoc_ProtTarMet$FDR < 0.05)
#	12,950/27,234 [47.6%]
sum(Assoc_ProtUntarMet$FDR < 0.05)
#	0/50,162 [0%]
sum(Assoc_ProtProt$Prot_FDR < 0.05)
#	50,398/71,022 [71%]
sum(Assoc_TarMetTarMet$TarMet_FDR < 0.05)
#	10,282/10,302 [99.8%]
sum(Assoc_UntarMetUntarMet$UntarMet_FDR < 0.05)
#	22,044/35,156 [62.7%]
sum(Assoc_TarMetUntarMet$UntarMet_FDR < 0.05)
#	0/19,176 [0%]

##	TOTAL (final) NUMBER OF ASSOCIATIONS|EDGES
#-------------------------------------------------------------------------
#		97,807,420 Total	out of	116,512,166	  [84%]
#		96,820,628	MolPhe1~MolPhe1	out of	113,398,462 	[85%]
#		986,792		MolPhe1~MolPhe2	out of	6,005,547   [16%]
#--------------------------------------------------------------------------------
##	SUBSET PAIRS SIGNIFICANT FOR ASSOCIATION
Assoc_ExprsProtAdj_sub <- subset(Assoc_ExprsProtAdj, subset=(FDR <= 0.05))
Assoc_ExprsTarMetAdj_sub <- subset(Assoc_ExprsTarMetAdj, subset=(FDR <= 0.05))
Assoc_ExprsUntarMetAdj_sub <- subset(Assoc_ExprsUntarMetAdj, subset=(FDR <= 0.05))
Assoc_ProtTarMetAdj_sub <- subset(Assoc_ProtTarMetAdj, subset=(FDR <= 0.05))
Assoc_ProtUntarMetAdj_sub <- subset(Assoc_ProtUntarMetAdj, subset=(FDR <= 0.05))
Assoc_ExprsExprsAdj_sub <- subset(Assoc_ExprsExprsAdj, subset=(Exprs_FDR <= 0.05))
Assoc_ProtProtAdj_sub <- subset(Assoc_ProtProtAdj, subset=(Prot_FDR <= 0.05))
Assoc_TarMetTarMetAdj_sub <- subset(Assoc_TarMetTarMetAdj, subset=(TarMet_FDR <= 0.05))
Assoc_UntarMetUntarMetAdj_sub <- subset(Assoc_UntarMetUntarMetAdj, subset=(UntarMet_FDR <= 0.05))
Assoc_TarMetUntarMetAdj_sub <- subset(Assoc_TarMetUntarMetAdj, subset=(UntarMet_FDR <= 0.05))

rownames(Assoc_ExprsProtAdj_sub) <- 1:nrow(Assoc_ExprsProtAdj_sub)
rownames(Assoc_ExprsTarMetAdj_sub) <- 1:nrow(Assoc_ExprsTarMetAdj_sub)
rownames(Assoc_ExprsUntarMetAdj_sub) <- 1:nrow(Assoc_ExprsUntarMetAdj_sub)
rownames(Assoc_ProtTarMetAdj_sub) <- 1:nrow(Assoc_ProtTarMetAdj_sub)
rownames(Assoc_ExprsExprsAdj_sub) <- 1:nrow(Assoc_ExprsExprsAdj_sub)
rownames(Assoc_ProtProtAdj_sub) <- 1:nrow(Assoc_ProtProtAdj_sub)
rownames(Assoc_TarMetTarMetAdj_sub) <- 1:nrow(Assoc_TarMetTarMetAdj_sub)
rownames(Assoc_UntarMetUntarMetAdj_sub) <- 1:nrow(Assoc_UntarMetUntarMetAdj_sub)

#---------------------------------------------------------
# CREATE SINGE FILE OF ASSOCIATIONS
#---------------------------------------------------------
Assoc_ExprsExprsAdj_sub$GeneNameA <- paste0(Assoc_ExprsExprsAdj_sub$GeneName1, "e")
Assoc_ExprsExprsAdj_sub$GeneNameB <- paste0(Assoc_ExprsExprsAdj_sub$GeneName2, "e")
Assoc_ExprsProtAdj_sub$GeneNameA <- paste0(Assoc_ExprsProtAdj_sub$GeneName, "e")
Assoc_ExprsProtAdj_sub$GeneNameB <- paste0(Assoc_ExprsProtAdj_sub$ProteinName, "p")
Assoc_ExprsTarMetAdj_sub$GeneNameA <- paste0(Assoc_ExprsTarMetAdj_sub$GeneName, "e")
Assoc_ExprsUntarMetAdj_sub$GeneNameA <- paste0(Assoc_ExprsUntarMetAdj_sub$GeneName, "e")
Assoc_ProtProtAdj_sub$GeneNameA <- paste0(Assoc_ProtProtAdj_sub$GeneName1, "p")
Assoc_ProtProtAdj_sub$GeneNameB <- paste0(Assoc_ProtProtAdj_sub$GeneName2, "p")
Assoc_ProtTarMetAdj_sub$GeneName2 <- paste0(Assoc_ProtTarMetAdj_sub$GeneName, "p")

Associations_ALL <- Assoc_ExprsProtAdj_sub[,c(10:11,5:8)]
colnames(Associations_ALL) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_ExprsExprs <- Assoc_ExprsExprsAdj_sub[,c(9:10,5:8)]
colnames(TEMP_Assoc_ExprsExprs) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_ProtProt <- Assoc_ProtProtAdj_sub[,c(10:11,5:8)]
colnames(TEMP_Assoc_ProtProt) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_ExprsTarMet <- Assoc_ExprsTarMetAdj_sub[,c(9,3:7)]
colnames(TEMP_Assoc_ExprsTarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_ExprsUntarMet <- Assoc_ExprsUntarMetAdj_sub[,c(10,3,5:8)]
colnames(TEMP_Assoc_ExprsUntarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_ProtTarMet <- Assoc_ProtTarMetAdj_sub[,c(9,3:7)]
colnames(TEMP_Assoc_ProtTarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_TarMetTarMet <- Assoc_TarMetTarMetAdj_sub[,c(1:6)]
colnames(TEMP_Assoc_TarMetTarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
TEMP_Assoc_UntarMetUntarMet <- Assoc_UntarMetUntarMetAdj_sub[,c(1,3,5:8)]
colnames(TEMP_Assoc_UntarMetUntarMet) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
Associations_ALL <- rbind(Associations_ALL, TEMP_Assoc_ExprsExprs, TEMP_Assoc_ProtProt, TEMP_Assoc_ExprsTarMet, TEMP_Assoc_ExprsUntarMet, TEMP_Assoc_ProtTarMet, TEMP_Assoc_TarMetTarMet, TEMP_Assoc_UntarMetUntarMet)

rownames(Associations_ALL) <- 1:nrow(Associations_ALL)
Associations_ALL$BETA <- round(Associations_ALL$BETA, digits=3)
Associations_ALL$T_value <- as.numeric(Associations_ALL$T_value)
Associations_ALL$T_value <- round(Associations_ALL$T_value, digits=3)
Associations_ALL$P_value <- round(Associations_ALL$P_value, digits=3)
Associations_ALL$FDR <- round(Associations_ALL$FDR, digits=3)
write.table(Associations_ALL, file = "ALL_ASSOC.tab", quote = FALSE, sep = "\t", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
Associations_ALL <- fread("ALL_ASSOC.tab")

### Omit Exprs~Exprs for network construction
tmp <- Assoc_ExprsProtAdj_sub[,c(10:11,5:8)]
colnames(tmp) <- c("MolPhe1","MolPhe2","BETA","T_value","P_value","FDR")
Associations_ALL2 <- rbind(tmp, TEMP_Assoc_ExprsTarMet, TEMP_Assoc_ExprsUntarMet, TEMP_Assoc_ProtProt, TEMP_Assoc_ProtTarMet, TEMP_Assoc_TarMetTarMet, TEMP_Assoc_UntarMetUntarMet)
rownames(Associations_ALL2) <- 1:nrow(Associations_ALL2)
Associations_ALL2$BETA <- round(Associations_ALL2$BETA, digits=3)
Associations_ALL2$T_value <- round(Associations_ALL2$T_value, digits=3)
Associations_ALL2$P_value <- round(Associations_ALL2$P_value, digits=3)
Associations_ALL2$FDR <- round(Associations_ALL2$FDR, digits=3)
write.table(Associations_ALL2, file = "ALL_ASSOC_minus.tab", quote = FALSE, sep = "\t", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
Associations_ALL2 <- fread("ALL_ASSOC_minus.tab")


#		TOTAL NUMBER OF POSSIBLE MOLECULES|NODES
length(unique(c(Associations_ALL$MolPhe1,Associations_ALL$MolPhe2))) #11,169
length(unique(c(Associations_ALL2$MolPhe1,Associations_ALL2$MolPhe2))) #11,166
#		10,643 + 267 + 103 + 188 === 11,201
#		11,166	/	11,201  [99.7%]


#		CREATE LABEL FILE
Assoc_Nodes_Labels <- data.frame(unique(c(Associations_ALL2$MolPhe1,Associations_ALL2$MolPhe2)))
colnames(Assoc_Nodes_Labels) <- "MolPhe"
#	11,166

lm_genes <- read.table( "../../2.ANALYSIS/2.Molecular_Phenotypes/3.GENES_DIRECT/lmResultsHbA1cAdjGenes_Age.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
lm_proteins <- read.table( "../../2.ANALYSIS/2.Molecular_Phenotypes/4.PROTEOMICS_DIRECT/lmResultsHbA1cAdjProteins_Age.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
lm_metabT <- read.table("../../2.ANALYSIS/2.Molecular_Phenotypes/5.METABOLITES_DIRECT/lmResultsHbA1cAdjTarMet_Age.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
lm_metabU <- fread("../../2.ANALYSIS/2.Molecular_Phenotypes/5.METABOLITES_DIRECT/lmResultsAdjUntarMet_HbA1cm.txt")
lm_metabU <- subset(lm_metabU, Age_Q<0.05)
#-------------------------------------------------------------------------------------
genes <- paste0(trimws(lm_genes$GeneName), "e")
proteins <- paste0(trimws(lm_proteins$GeneName), "p")
metabolitesT <- lm_metabT$Metabolite
metabolitesU <- lm_metabU$BIOCHEMICAL2
#-------------------------------------------------------------------------------------------------------------------
for(i in 1:nrow(Assoc_Nodes_Labels)){
  if(Assoc_Nodes_Labels$MolPhe[i] %in% genes ){
    Assoc_Nodes_Labels$Name[i] <- as.character(strsplit(Assoc_Nodes_Labels$MolPhe[i], "e"))
    Assoc_Nodes_Labels$Label[i] <- "Expression"}
  else if(Assoc_Nodes_Labels$MolPhe[i] %in% proteins){
    Assoc_Nodes_Labels$Name[i] <- as.character(strsplit(Assoc_Nodes_Labels$MolPhe[i], "p"))
    Assoc_Nodes_Labels$Label[i] <- "Protein"}
  else if(Assoc_Nodes_Labels$MolPhe[i] %in% metabolitesT){Assoc_Nodes_Labels$Name[i] <- with(lm_metabT, CommonName[Metabolite==Assoc_Nodes_Labels$MolPhe[i]])
  Assoc_Nodes_Labels$Label[i] <- "Metabolite"}
  else if(Assoc_Nodes_Labels$MolPhe[i] %in% metabolitesU){
    Assoc_Nodes_Labels$Name[i] <- Assoc_Nodes_Labels$MolPhe[i]
    Assoc_Nodes_Labels$Label[i] <- "Metabolite"}}
write.table(Assoc_Nodes_Labels, file = "ALL_NODES_Labels.tab", quote = FALSE, sep = "\t", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
Assoc_Nodes_Labels <- fread("ALL_NODES_Labels.tab")


#		INCLUDE BIC SCORE FOR EACH PAIR (BAYES)
BIC_Exprs_Exprs <- fread("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/9.RELATIONS/RESULTS_BAYES_AGE/BIC_Exprs_Exprs_Annotated.tab")
BIC_Exprs_Prot <- fread("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/9.RELATIONS/RESULTS_BAYES_AGE/BIC_Exprs_Prot_Annotated.tab")
BIC_Exprs_TarMet <- fread("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/9.RELATIONS/RESULTS_BAYES_AGE/BIC_Exprs_TarMet_Annotated.tab")
BIC_Exprs_UntarMet <- fread("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/9.RELATIONS/RESULTS_BAYES_AGE/BIC_Exprs_UntarMet_Annotated.tab")
BIC_Prot_Prot <- fread("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/9.RELATIONS/RESULTS_BAYES_AGE/BIC_Prot_Prot_Annotated.tab")
BIC_Prot_TarMet <- fread("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/9.RELATIONS/RESULTS_BAYES_AGE/BIC_Prot_TarMet_Annotated.tab")
BIC_TarMet_TarMet <- fread("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/9.RELATIONS/RESULTS_BAYES_AGE/BIC_TarMet_TarMet_Annotated.tab")
BIC_UntarMet_UntarMet <- fread("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/9.RELATIONS/RESULTS_BAYES_AGE/BIC_UntarMet_UntarMet_Annotated.tab")


#     SUBSET FILES FOR BIC >= 6
BIC_Exprs_Exprs <- subset(BIC_Exprs_Exprs, SIGNIF != "*")
BIC_Exprs_Prot <- subset(BIC_Exprs_Prot, SIGNIF != "*")
BIC_Exprs_TarMet <- subset(BIC_Exprs_TarMet, SIGNIF != "*")
BIC_Exprs_UntarMet <- subset(BIC_Exprs_UntarMet, SIGNIF != "*")
BIC_Prot_Prot <- subset(BIC_Prot_Prot, SIGNIF != "*")
BIC_Prot_TarMet <- subset(BIC_Prot_TarMet, SIGNIF != "*")
BIC_TarMet_TarMet <- subset(BIC_TarMet_TarMet, SIGNIF != "*")
BIC_UntarMet_UntarMet <- subset(BIC_UntarMet_UntarMet, SIGNIF != "*")

# Total:88,153,775
sum(BIC_Exprs_Exprs$TOP == "model2")	#modelGENE: 2 -> 1
#[1]	8,820,297 [10%]
sum(BIC_Exprs_Exprs$TOP == "model1")	#modelPROTEIN: 1 -> 2
#[1]	9,998,129 [11%]
sum(BIC_Exprs_Exprs$TOP == "modelINDEP")
#[1]	69,335,349 [79%]

# Total: 460,694
sum(BIC_Exprs_Prot$TOP == "model2")	#modelGENE: 2 -> 1
#[1]	104,567 [23%]
sum(BIC_Exprs_Prot$TOP == "model1")	#modelPROTEIN: 1 -> 2
#[1]	303,892 [66%]
sum(BIC_Exprs_Prot$TOP == "modelINDEP")
#[1]	52,235 [11%]

# Total: 89,377
sum(BIC_Exprs_TarMet$TOP == "model2")	#modelGENE: 2 -> 1
#[1]	46,339 [52%]
sum(BIC_Exprs_TarMet$TOP == "model1")	#modelPROTEIN: 1 -> 2
#[1]	36,181 [40%]
sum(BIC_Exprs_TarMet$TOP == "modelINDEP")
#[1]	6,857 [8%]

# Total: 47,788
sum(BIC_Exprs_UntarMet$TOP == "model2")	#modelGENE: 2 -> 1
#[1]	19,014 [40%]
sum(BIC_Exprs_UntarMet$TOP == "model1")	#modelPROTEIN: 1 -> 2
#[1]	22,299 [47%]
sum(BIC_Exprs_UntarMet$TOP == "modelINDEP")
#[1]	6,475 [14%]

# Merged: Exprs ~ Met
# model2:       65,353/137,165 --> 48%
# model1:       58,480/137,165 --> 43%
# independent:  13,332/137,165 --> 8%

# Total: 62,530
sum(BIC_Prot_Prot$TOP == "model2")	#modelGENE: 2 -> 1
#[1]	8,063 [13%]
sum(BIC_Prot_Prot$TOP == "model1")	#modelPROTEIN: 1 -> 2
#[1]	10,948 [18%]
sum(BIC_Prot_Prot$TOP == "modelINDEP")
#[1]	43,519 [70%]

# Total: 8,943
sum(BIC_Prot_TarMet$TOP == "model2")	#modelGENE: 2 -> 1
#[1]	5,032 [56%]
sum(BIC_Prot_TarMet$TOP == "model1")	#modelPROTEIN: 1 -> 2
#[1]	2,140 [24%]
sum(BIC_Prot_TarMet$TOP == "modelINDEP")
#[1]	1,771 [20%]

# Total: 10,204
sum(BIC_TarMet_TarMet$TOP == "model2")	#modelGENE: 2 -> 1
#[1]	86 [0.8%]
sum(BIC_TarMet_TarMet$TOP == "model1")	#modelPROTEIN: 1 -> 2
#[1]	200 [2%]
sum(BIC_TarMet_TarMet$TOP == "modelINDEP")
#[1]	9,918 [97%]

# Total: 2,940
sum(BIC_UntarMet_UntarMet$TOP == "model2")	#modelGENE: 2 -> 1
#[1]	0
sum(BIC_UntarMet_UntarMet$TOP == "model1")	#modelPROTEIN: 1 -> 2
#[1]	20 [0.7%]
sum(BIC_UntarMet_UntarMet$TOP == "modelINDEP")
#[1]	2,920 [99%]

# Merged: Met ~ Met
# model2:       86/13,144 --> 0.7%
# model1:       220/13,144 --> 1.7%
# independent:  12,838/13,144 --> 98%

#Assoc_Nodes_Labels <- fread("ALL_NODES_Labels.tab") 
for(i in 1:nrow(BIC_Exprs_UntarMet)){BIC_Exprs_UntarMet$MolPhe2[i] <- with(Assoc_Nodes_Labels, MolPhe[Name==BIC_Exprs_UntarMet$MolPhe2[i]])}
for(i in 1:nrow(BIC_UntarMet_UntarMet)){
  BIC_UntarMet_UntarMet$MolPhe1[i] <- with(Assoc_Nodes_Labels, MolPhe[Name==BIC_UntarMet_UntarMet$MolPhe1[i]])
  BIC_UntarMet_UntarMet$MolPhe2[i] <- with(Assoc_Nodes_Labels, MolPhe[Name==BIC_UntarMet_UntarMet$MolPhe2[i]])}

# CREATE SINGE FILE OF BIC
#-------------------------------------------
BIC_ALL <- BIC_Exprs_Prot[,c(2,4:7)]
colnames(BIC_ALL) <- c("MolPhe1","MolPhe2","Direction","Weight","Weight_Cat")
TEMP_BIC_ExprsExprs <- BIC_Exprs_Exprs[,c(2,4:7)]
colnames(TEMP_BIC_ExprsExprs) <- c("MolPhe1","MolPhe2","Direction","Weight","Weight_Cat")
TEMP_BIC_ProtProt <- BIC_Prot_Prot[,c(2,4:7)]
colnames(TEMP_BIC_ProtProt) <- c("MolPhe1","MolPhe2","Direction","Weight","Weight_Cat")
TEMP_BIC_ExprsTarMet <- BIC_Exprs_TarMet[,c(2:3,5:7)]
colnames(TEMP_BIC_ExprsTarMet) <- c("MolPhe1","MolPhe2","Direction","Weight","Weight_Cat")
TEMP_BIC_ExprsUntarMet <- BIC_Exprs_UntarMet[,c(2:6)]
colnames(TEMP_BIC_ExprsUntarMet) <- c("MolPhe1","MolPhe2","Direction","Weight","Weight_Cat")
TEMP_BIC_ProtTarMet <- BIC_Prot_TarMet[,c(2:3,5:7)]
colnames(TEMP_BIC_ProtTarMet) <- c("MolPhe1","MolPhe2","Direction","Weight","Weight_Cat")
TEMP_BIC_TarMetTarMet <- BIC_TarMet_TarMet[,c(1,3,5:7)]
colnames(TEMP_BIC_TarMetTarMet) <- c("MolPhe1","MolPhe2","Direction","Weight","Weight_Cat")
TEMP_BIC_UntarMetUntarMet <- BIC_UntarMet_UntarMet
colnames(TEMP_BIC_UntarMetUntarMet) <- c("MolPhe1","MolPhe2","Direction","Weight","Weight_Cat")

BIC_ALL_T <- rbind(BIC_ALL, TEMP_BIC_ExprsExprs, TEMP_BIC_ExprsTarMet, TEMP_BIC_ExprsUntarMet, TEMP_BIC_ProtProt, TEMP_BIC_ProtTarMet, TEMP_BIC_TarMetTarMet, TEMP_BIC_UntarMetUntarMet)
# 88,836,251
rownames(BIC_ALL_T) <- 1:nrow(BIC_ALL_T)
BIC_ALL_T$Weight <- round(BIC_ALL_T$Weight, digits=3)
write.table(BIC_ALL_T, file = "ALL_BIC.tab", quote = FALSE, sep = "\t", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
BIC_ALL_T <- fread("S70_Relationships_AgeRelatedMolecules_BIC.tab")

BIC_ALL_M <- rbind(BIC_ALL, TEMP_BIC_ExprsTarMet, TEMP_BIC_ExprsUntarMet, TEMP_BIC_ProtProt, TEMP_BIC_ProtTarMet, TEMP_BIC_TarMetTarMet, TEMP_BIC_UntarMetUntarMet)
#	875,597
rownames(BIC_ALL_M) <- 1:nrow(BIC_ALL_M)
BIC_ALL_M$Weight <- round(BIC_ALL_M$Weight, digits=3)
write.table(BIC_ALL_M, file = "ALL_BIC_minus.tab", quote = FALSE, sep = "\t", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
BIC_ALL_M <- fread("ALL_BIC_minus.tab")


#     MAKE FILE WITH ASSOCIATIONS AND DIRECTION COLUMN
Associations_ALL2$Pair <- paste(Associations_ALL2$MolPhe1,Associations_ALL2$MolPhe2, sep="_")
BIC_ALL_M$Pair <- paste(BIC_ALL_M$MolPhe1,BIC_ALL_M$MolPhe2, sep="_")

# 1,069,516-682,476 = 387,040 (no Bayes info)
length(intersect(BIC_ALL_M$Pair, Associations_ALL2$Pair))
# [1] 513,067
# 682,476 - 513,067 = 169,409
index <- which(BIC_ALL_M$Pair %!in% Associations_ALL2$Pair)
length(index)
#[1] 168,261
which(Associations_ALL2$Pair == "RP11-206L10.2e_REG1Ap")
Associations_ALL2[2538,]

for(i in 1650:nrow(Associations_ALL2)){
  if(Associations_ALL2$Pair[i] %in% BIC_ALL_M$Pair){
    Associations_ALL2$Direction[i] <- with(BIC_ALL_M, Direction[Pair == Associations_ALL2$Pair[i]])
    Associations_ALL2$Weight[i] <- with(BIC_ALL_M, Weight[Pair == Associations_ALL2$Pair[i]])
    Associations_ALL2$WeightCAT[i] <- with(BIC_ALL_M, Weight_Cat[Pair == Associations_ALL2$Pair[i]])}
  else{
    print(i)
    Associations_ALL2$Direction[i] <- NA
    Associations_ALL2$Weight[i] <-NA
    Associations_ALL2$WeightCAT[i] <- NA}}

write.table(Associations_ALL2, file = "S78_Relationships_AgeRelatedMolecules_BIC.tab", quote = FALSE, sep = "\t", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
network <- fread("S78_Relationships_AgeRelatedMolecules_BIC.tab")

##########################################################################################################################################################################
##########################################################################################################################################################################
#                                           CYTOSCAPE - mTOR
##########################################################################################################################################################################
# 5,071 of 11,166 nodes [45%] connected with >100 nodes (degree)
# Top degree: The most central element of this age-network was protein IL6 connected with 9,015 molecules, followed by proteins LILRA5 and PLAUR (degree > 8000).
##########################################################################################################################################################################
# 1. Import "S70_NETWORK_DATA.tabb"
# 2. Import "S70_NODES_Labels" [table from file]
# 3. Subset based on mTOR genes and proteins [upper right coner with spaces in between]   --	35/104 in pathway
# 4. Select 1st neighbours of selected & edges between them -> new network: 759 nodes, 112,758 edges
# 5. (+) Add "Column Filter" - Adjust limit of FDR @ 1.0E-10
# 6. (+) Add "Column Filter" - Remove edges with NA OR modelINDEP in Direction ---> 647 nodes, 13,899 edges
# -	Layout > Degree sorted circle layout - selected nodes only (mTOR)   ---  move circle
# 7. Export node table: mTOR_nodes.csv
--------------------------------------------------------------------------------------------------------------------------------------
#  WE NEED TO ADD BETA SIGN COLUMN FOR SEX [R]
lm_genes <- read.delim("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/3.GENES_DIRECT/lmResultsAdjGenes_HbA1cm.txt", header = TRUE)
lm_proteins <- read.delim("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/4.PROTEOMICS_DIRECT/lmResultsAdjProteins_HbA1cm.txt", header = TRUE)
lm_metabT <- read.delim("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/5.METABOLITES_DIRECT/lmResultsAdjTarMet_HbA1cm.txt", header = TRUE)
lm_metabU <- read.delim("~/Desktop/PhD/2.ANALYSIS/2.Molecular_Phenotypes/5.METABOLITES_DIRECT/lmResultsAdjUntarMet_HbA1cm.txt", header = TRUE)
lm_genes$GeneName <- trimws(lm_genes$GeneName)
lm_proteins$GeneName <- trimws(lm_proteins$GeneName)
lm_metabT$CommonName <- trimws(lm_metabT$CommonName)
lm_metabU$BIOCHEMICAL2 <- trimws(lm_metabU$BIOCHEMICAL2)
network <- read.table("mTOR_nodes.csv",header = TRUE, sep=",")

sex_mol <- c(paste0(subset(lm_genes,select=GeneName,subset=Sex_Q < 0.05)[,1],"e"),paste0(subset(lm_proteins,select=GeneName,subset=Sex_Q < 0.05)[,1],"p"),subset(lm_metabT,select=Metabolite ,subset=Sex_Q < 0.05)[,1],subset(lm_metabU,select=BIOCHEMICAL2,subset=Sex_Q < 0.05)[,1])
for(i in 1:nrow(network)){
  if(network$shared.name[i] %in% sex_mol){network$Sex[i] <-  "yes"}
  else{network$Sex[i] <- "no"}}
for(i in 1:nrow(network)){
  print(i)
  if(network$Label[i] == "Expression"){network$SexBeta[i] <-  with(lm_genes, SexBetas[GeneName == network$name[i]])}
  else	if(network$Label[i] == "Protein"){network$SexBeta[i] <-  with(lm_proteins, SexBetas[GeneName == network$name[i]])}
  else	if(network$Label[i] == "Metabolite"){
    if(network$shared.name[i] %in% lm_metabT$Metabolite){network$SexBeta[i] <-  with(lm_metabT, SexBetas[CommonName == network$name[i]])}
    else if(network$shared.name[i] %in% lm_metabU$BIOCHEMICAL2){network$SexBeta[i] <-  with(lm_metabU, SexBetas[BIOCHEMICAL2 == network$name[i]])}}}
network$SexBeta <- trimws(network$SexBeta)

network$Label_Sex_Beta <- paste(network$Label,network$Sex,network$SexBeta,sep="_")

mTOR <- c("AKT1e","AKT2e","AKT3e","BRAFe","CAB39e","CAB39Le","DDIT4e","EIF4Be","EIF4Ee","EIF4E1Be","EIF4E2e","EIF4EBP1e","HIF1Ae",
          "IGF1e","INSe","MAPK1e","MAPK3e","MLST8e","MTORe","PDPK1e","PGFe","PIK3CAe","PIK3CBe","PIK3CDe","PIK3CGe","PIK3R1e",
          "PIK3R2e","PIK3R3e","PIK3R5e","PRKAA1e","PRKAA2e","RHEBe","RICTORe","RPS6e","RPS6KA1e","RPS6KA2e","RPS6KA3e",
          "RPS6KA6e","RPS6KB1e","RPS6KB2e","RPTORe","STK11e","STRADAe","TSC1e","TSC2e","ULK1e","ULK2e","ULK3e","VEGFAe","VEGFBe",
          "VEGFCe","VEGFDe","AKT1p","AKT2p","AKT3p","BRAFp","CAB39p","CAB39Lp","DDIT4p","EIF4Bp","EIF4Ep","EIF4E1Bp","EIF4E2p",
          "EIF4EBP1p","HIF1Ap","IGF1p","INSp","MAPK1p","MAPK3p","MLST8p","MTORp","PDPK1p","PGFp","PIK3CAp","PIK3CBp","PIK3CDp",
          "PIK3CGp","PIK3R1p","PIK3R2p","PIK3R3p","PIK3R5p","PRKAA1p","PRKAA2p","RHEBp","RICTORp","RPS6p","RPS6KA1p","RPS6KA2p",
          "RPS6KA3p","RPS6KA6p","RPS6KB1p","RPS6KB2p","RPTORp","STK11p","STRADAp","TSC1p","TSC2p","ULK1p","ULK2p",
          "ULK3p","VEGFAp","VEGFBp","VEGFCp","VEGFDp")
for(i in 1:nrow(network)){
  if(network$shared.name[i] %in% mTOR){network$mTOR[i] <-  "yes"}
  else{network$mTOR[i] <- "no"}}
# REMOVE UNNECESSARY/DUPLICATED COLUMNS!!! 
write.table(network, file = "mTOR_nodes_sexBeta.txt", quote = FALSE, sep = "\t", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
----------------------------------------------------------------------------------------------------------------------------------------
# 7. Import "mTOR_nodes_sexBeta.txt" [table from file - Import Data as Node Table Columns]
# 8. Use "Label_Sex_Beta" to colour nodes
# 9. Styles > Edges > Source Arrow Shape (Direction): model2 = Arrow
#                     Target Arrow Shape (Direction): model1 = Arrow
# 10.Select 35 mTOR + nodes connected to edges
# 11. Import to new network --->	237 nodes, 557 edges
# 12. Tools > Analyse Network (undirected)
# 13. Remove UNKNOWN metabolites
# 14. Remove Proteins with degree<=3  
# 15. Remove Metabolites with degree<=2 
# 16. Tweak visuals    
# 17. Recalculate degree
# 17. Export & save as .cys
##########################################################################################################################################################################
##########################################################################################################################################################################
#                                           CYTOSCAPE - T2D
##########################################################################################################################################################################
##########################################################################################################################################################################
# 1. Import "S70_NETWORK_DATA.tab"
# 2. Import NEW Network (ALL_SNPs)
# 3. Tools > Merge networks > Union		BY SHARED NAME!
# 4. Import "S70_NODES_SNPs_Labels.tab" [table from file]
# 5. Subset based on  T2D coloc hits  --	7/12 (unique) 
# 6. Select 1st neighbours of selected & edges between them [x2 to get SNPs] -> new network: 350 nodes, 43,630 edges
# 7. (+) Add "Column Filter" - Adjust limit of FDR @ 0.001
# 8. (+) Add "Column Filter" - Remove edges with NA OR modelINDEP in Direction 
# 9. Select nodes connected by selected edges
# 10. Manually select remaining SNPs and their edges!
# 11. Export to new network: 299 nodes, 9248 edges
# 12. Recalculate degree (Tools > Analyse network)
# 14. ALL PROTEINS HAVE DEGREE 17 - 270 -- Remove Proteins with degree<90
# 15. ALL METABOLITES HAVE DEGREE 1-74 -- Remove Metabolites with degree<65		------------> 67 nodes, 1,211 edges

