# THIS FILE CONTAINS THE R CODE USED TO IDENTIFY RELATIONSHIPS (ASSOCIATIONS) BETWEEN
# AGE-RELATED MOLECULAR PHENOTYPES OF THE SAME TYPE [GENES, PROTEINS AND METABOLITES (TARGETED & UNTARGETED DATASET)]
# USING LINEAR REGRESSION MODELS
# FOR ASSOCIATIONS BETWEEN GENE EXPRESSION MOLECULES, AN EXTERNAL BASH SCRIPT WAS USED
# MULTIPLE CORRECTION METHODS WAS APPLIED FOR THE SUBSET OF PAIRS WHERE BOTH MOLECULES WERE ASSOCIATED WITH AGE
# (FOR SUBSEQUENT CAUSAL INFERENCE ANALYSIS)
# BAYESIAN CAUSAL INFERENCE WAS PERFORMED USING BASH SCRIPTS
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#		MOLECULAR PHENOTYPE RELATIONSHIPS (ASSOCIATIONS): LINEAR REGRESSION - GET RESULTS [MolPheno1 ~	 MolPheno1]
##	RUN		EXPS ~ EXPRS

index_age_genes <- which(as.data.frame(resid_genes_TechBiol)[,2] %in% res_GenesHbA1cAdj_Age$GeneName)
resid_genes_TechBiol_age <- subset(resid_genes_TechBiol, subset=(GeneID %in% res_GenesHbA1cAdj_Age$GeneID))
write.table(resid_genes_TechBiol_age, file = "./RELATIONS/ExprsResiduals_AgeSUBSET.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(covar_Exprsf$Age, file = "./RELATIONS/Exprs_Age3027.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

#-------------------------------------------------------------------------------------
# qsub  /home/dafmic/CrossSectional/RELATIONS/ExprsExprs.sh
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#echo "cat /users/home/dafmic/CrossSectional/RELATIONS/results2/ExprsExprs_Assoc_{1..10643} | gzip -c > /users/home/dafmic/CrossSectional/RELATIONS/results2/ExprsExprs_Assoc_ALL.txt.gz"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# IMPORT		EXPRS ~ EXPRS
relat_ExprsExprs <- read.table("/users/home/dafmic/CrossSectional/RELATIONS/results2/ExprsExprs_Assoc_ALL.txt.gz", header = FALSE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
is.data.frame(relat_ExprsExprs)
colnames(relat_ExprsExprs) <- c("GeneID1","GeneName1","GeneID2","GeneName2","Exprs","Exprs_Tvalue","Exprs_pvalue")

# RUN		PROT ~ PROT
relat_ProtProt <- lmMolPheAGEx2(exprs1 = resid_proteins_TechBiol[,3:3025],
																	exprs2 = resid_proteins_TechBiol[,3:3025],
																	exprs1_annot = resid_proteins_TechBiol[,1:2],
																	exprs2_annot = resid_proteins_TechBiol[,1:2],
																	covar = covar_Protf$Age,
																	exprs2Name = "Prot")															
relat_ProtProt <- as.data.frame(relat_ProtProt)
colnames(relat_ProtProt)[1:4] <- c("ProteinID1","GeneName1","ProteinID2","GeneName2")
# RUN		MET_TAR ~ MET_TAR
relat_TarMetTarMet <- lmMolPheAGEx2(exprs1 = resid_tarMet_TechBiol[,2:3028],
																	exprs2 = resid_tarMet_TechBiol[,2:3028],
																	exprs1_annot = as.data.frame(resid_tarMet_TechBiol[,1]),
																	exprs2_annot = as.data.frame(resid_tarMet_TechBiol[,1]),
																	covar = covar_Targf$Age,
																	exprs2Name = "TarMet")															
relat_TarMetTarMet <- as.data.frame(relat_TarMetTarMet)
colnames(relat_TarMetTarMet)[1:2] <- c("MetaboliteID1","MetaboliteID2")
# RUN		MET_UNTAR ~ MET_UNTAR
relat_UntarMetUntarMet <- lmMolPheAGEx2(exprs1 = resid_untarMet_TechBiol[,3:2998],
																	exprs2 = resid_untarMet_TechBiol[,3:2998],
																	exprs1_annot = resid_untarMet_TechBiol[,1:2],
																	exprs2_annot = resid_untarMet_TechBiol[,1:2],
																	covar = covar_Untargf$Age,
																	exprs2Name = "UntarMet")															
relat_UntarMetUntarMet <- as.data.frame(relat_UntarMetUntarMet)
colnames(relat_UntarMetUntarMet)[1:4] <- c("MetaboliteID1","SuperPathway1","MetaboliteID2","SuperPathway2")
# RUN		MET_TAR ~ MET_UNTAR
relat_TarMetUntarMet <- lmMolPheAGEx2(exprs1 = resid_tarMet_TechBiol_sub3023[,2:2997],
																		exprs2 = resid_untarMet_TechBiol[,3:2998],
																		exprs1_annot = as.data.frame(resid_tarMet_TechBiol_sub3023[,1]),
																		exprs2_annot = resid_untarMet_TechBiol[,1:2],
																		covar = covar_Untargf$Age,
																		exprs2Name = "UntarMet")															
relat_TarMetUntarMet <- as.data.frame(relat_TarMetUntarMet)
colnames(relat_TarMetUntarMet)[1:3] <- c("MetaboliteID1","MetaboliteID2","SuperPathway2")


#		SUBSET SIGNIFICANT PAIRS OF AGE-RELATED MOLECULES: MolPheno1 ~ MolPheno2
res_GenesHbA1cAdj_Age <- read.table( "/users/home/dafmic/CrossSectional/Transcriptomics/lmResultsAdjGenes_HbA1cm.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
res_ProteinsHbA1cAdj_Age <- read.table( "/users/home/dafmic/CrossSectional/Proteomics/lmResultsAdjProteins_HbA1cm.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
res_TargHbA1cAdj_Age <- read.table("/users/home/dafmic/CrossSectional/Metabolomics/lmResultsAdjTarMet_HbA1cm.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
res_UntargHbA1cAdj_Age <- read.table("/users/home/dafmic/CrossSectional/Metabolomics/lmResultsAdjUntarMet_HbA1cm.txt",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)

relat_ExprsProt$GeneID <- trimws(relat_ExprsProt$GeneID)
relat_ExprsTarMet$GeneID <- trimws(relat_ExprsTarMet$GeneID)
relat_ExprsUntarMet$GeneID <- trimws(relat_ExprsUntarMet$GeneID)
relat_ProtTarMet$ProteinID <- trimws(relat_ProtTarMet$ProteinID)
relat_ProtUntarMet$ProteinID <- trimws(relat_ProtUntarMet$ProteinID)
relat_ExprsProt$ProteinID <- trimws(relat_ExprsProt$ProteinID)
relat_ExprsTarMet$MetaboliteID <- trimws(relat_ExprsTarMet$MetaboliteID)
relat_ExprsUntarMet$BIOCHEMICAL2 <- trimws(relat_ExprsUntarMet$BIOCHEMICAL2)
relat_ProtTarMet$MetaboliteID <- trimws(relat_ProtTarMet$MetaboliteID)
relat_ProtUntarMet$BIOCHEMICAL2 <- trimws(relat_ProtUntarMet$BIOCHEMICAL2)
res_GenesHbA1cAdj_Age$GeneID <- trimws(res_GenesHbA1cAdj_Age$GeneID)
res_ProteinsHbA1cAdj_Age$ProteinID <- trimws(res_ProteinsHbA1cAdj_Age$ProteinID)
res_TargHbA1cAdj_Age$Metabolite <- trimws(res_TargHbA1cAdj_Age$Metabolite)
res_UntargHbA1cAdj_Age$BIOCHEMICAL2 <- trimws(res_UntargHbA1cAdj_Age$BIOCHEMICAL2)

##	SUBSET ASSOCIATIONS WHERE BOTH MOLPHENO ARE ASSOC. WITH AGE
relat_ExprsProtAdj <- subset(relat_ExprsProt, subset=((GeneID %in% res_GenesHbA1cAdj_Age$GeneID)&(ProteinID %in% res_ProteinsHbA1cAdj_Age$ProteinID)))
relat_ExprsTarMetAdj <- subset(relat_ExprsTarMet, subset=((GeneID %in% res_GenesHbA1cAdj_Age$GeneID)&(MetaboliteID %in% res_TargHbA1cAdj_Age$Metabolite)))
relat_ExprsUntarMetAdj <- subset(relat_ExprsUntarMet, subset=((GeneID %in% res_GenesHbA1cAdj_Age$GeneID)&(BIOCHEMICAL2 %in% res_UntargHbA1cAdj_Age$BIOCHEMICAL2)))
relat_ProtTarMetAdj <- subset(relat_ProtTarMet, subset=((ProteinID %in% res_ProteinsHbA1cAdj_Age$ProteinID)&(MetaboliteID %in% res_TargHbA1cAdj_Age$Metabolite)))
relat_ProtUntarMetAdj <- subset(relat_ProtUntarMet, subset=((ProteinID %in% res_ProteinsHbA1cAdj_Age$ProteinID)&(BIOCHEMICAL2 %in% res_UntargHbA1cAdj_Age$BIOCHEMICAL2)))

##	MULTIPLE CORRECTION
relat_ExprsProtAdj$Prot_Q <- qvalue(as.numeric(relat_ExprsProtAdj[,"Prot_pvalue"]))$qvalues
relat_ExprsTarMetAdj$TarMet_Q <- qvalue(as.numeric(relat_ExprsTarMetAdj[,"TarMet_pvalue"]))$qvalues
relat_ExprsUntarMetAdj$UntarMet_Q <- qvalue(as.numeric(relat_ExprsUntarMetAdj[,"UntarMet_pvalue"]))$qvalues
relat_ProtTarMetAdj$TarMet_Q <- qvalue(as.numeric(relat_ProtTarMetAdj[,"TarMet_pvalue"]))$qvalues
relat_ProtUntarMetAdj$UntarMet_Q <- qvalue(as.numeric(relat_ProtUntarMetAdj[,"UntarMet_pvalue"]))$qvalues
relat_ExprsProtAdj$Prot_FDR <- p.adjust(as.numeric(relat_ExprsProtAdj$Prot_pvalue), method = "BH")
relat_ExprsTarMetAdj$TarMet_FDR <- p.adjust(as.numeric(relat_ExprsTarMetAdj$TarMet_pvalue), method = "BH")
relat_ExprsUntarMetAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_ExprsUntarMetAdj$UntarMet_pvalue), method = "BH")
relat_ProtTarMetAdj$TarMet_FDR <- p.adjust(as.numeric(relat_ProtTarMetAdj$TarMet_pvalue), method = "BH")
relat_ProtUntarMetAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_ProtUntarMetAdj$UntarMet_pvalue), method = "BH")

write.table(relat_ExprsProtAdj, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ExprsProt_Assoc.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_ExprsTarMetAdj, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ExprsTarMet_Assoc.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_ExprsUntarMetAdj, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ExprsUntarMet_Assoc.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_ProtTarMetAdj, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ProtTarMet_Assoc.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_ProtUntarMetAdj, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ProtUntarMet_Assoc.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

sum(relat_ExprsProtAdj$Prot_Q < 0.05)
#	1,248,850
sum(relat_ExprsTarMetAdj$TarMet_Q < 0.05)
#	281,670
sum(relat_ExprsUntarMetAdj$UntarMet_Q < 0.05)
#	136,148
sum(relat_ProtTarMetAdj$TarMet_Q < 0.05)
#	25,257
sum(relat_ProtUntarMetAdj$UntarMet_Q < 0.05)
#	0
sum(relat_ExprsProtAdj$Prot_FDR < 0.05)
#	1,036,750
sum(relat_ExprsTarMetAdj$TarMet_FDR < 0.05)
#	217,969
sum(relat_ExprsUntarMetAdj$UntarMet_FDR < 0.05)
#	111,943
sum(relat_ProtTarMetAdj$TarMet_FDR < 0.05)
#	19,852
sum(relat_ProtUntarMetAdj$UntarMet_FDR < 0.05)
#	0

##	SUBSET PAIRS SIGNIFICANT FOR ASSOCIATION
relat_ExprsProtAdj_sub <- subset(relat_ExprsProtAdj, subset=(Prot_FDR <= 0.05))
#	1,036,750
relat_ExprsTarMetAdj_sub <- subset(relat_ExprsTarMetAdj, subset=(TarMet_FDR <= 0.05))
#	217,969
relat_ExprsUntarMetAdj_sub <- subset(relat_ExprsUntarMetAdj, subset=(UntarMet_FDR <= 0.05))
#	111,943
relat_ProtTarMetAdj_sub <- subset(relat_ProtTarMetAdj, subset=(TarMet_FDR <= 0.05))
#	19,852
relat_ProtUntarMetAdj_sub <- subset(relat_ProtUntarMetAdj, subset=(UntarMet_FDR <= 0.05))
#	0

rownames(relat_ExprsProtAdj_sub) <- 1:nrow(relat_ExprsProtAdj_sub)
rownames(relat_ExprsTarMetAdj_sub) <- 1:nrow(relat_ExprsTarMetAdj_sub)
rownames(relat_ExprsUntarMetAdj_sub) <- 1:nrow(relat_ExprsUntarMetAdj_sub)
rownames(relat_ProtTarMetAdj_sub) <- 1:nrow(relat_ProtTarMetAdj_sub)
#rownames(relat_ProtUntarMetAdj_sub) <- 1:nrow(relat_ProtUntarMetAdj_sub)

write.table(relat_ExprsProtAdj_sub, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ExprsProt_Assoc_Signif.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_ExprsTarMetAdj_sub, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ExprsTarMet_Assoc_Signif.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_ExprsUntarMetAdj_sub, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ExprsUntarMet_Assoc_Signif.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_ProtTarMetAdj_sub, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ProtTarMet_Assoc_Signif.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
#write.table(relat_ProtUntarMetAdj_sub, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ProtUntarMet_Assoc_Signif.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)


#		SUBSET SIGNIFICANT PAIRS OF AGE-RELATED MOLECULES: MolPheno1 ~ MolPheno1
relat_ExprsExprs$GeneID2 <- trimws(relat_ExprsExprs$GeneID2)
relat_ProtProt$ProteinID1 <- trimws(relat_ProtProt$ProteinID1)
relat_ProtProt$ProteinID2 <- trimws(relat_ProtProt$ProteinID2)
relat_TarMetTarMet$MetaboliteID1 <- trimws(relat_TarMetTarMet$MetaboliteID1)
relat_TarMetTarMet$MetaboliteID2 <- trimws(relat_TarMetTarMet$MetaboliteID2)
relat_UntarMetUntarMet$MetaboliteID1 <- trimws(relat_UntarMetUntarMet$MetaboliteID1)
relat_UntarMetUntarMet$MetaboliteID2 <- trimws(relat_UntarMetUntarMet$MetaboliteID2)
relat_TarMetUntarMet$MetaboliteID1 <- trimws(relat_TarMetUntarMet$MetaboliteID1)
relat_TarMetUntarMet$MetaboliteID2 <- trimws(relat_TarMetUntarMet$MetaboliteID2)

relat_ExprsExprsAdj <- subset(relat_ExprsExprs, subset=((GeneID1 %in% res_GenesHbA1cAdj_Age$GeneID)&(GeneID2 %in% res_GenesHbA1cAdj_Age$GeneID)))
relat_ProtProtAdj <- subset(relat_ProtProt, subset=((ProteinID1 %in% res_ProteinsHbA1cAdj_Age$ProteinID)&(ProteinID2 %in% res_ProteinsHbA1cAdj_Age$ProteinID)))
relat_TarMetTarMetAdj <- subset(relat_TarMetTarMet, subset=((MetaboliteID1 %in% res_TargHbA1cAdj_Age$Metabolite)&(MetaboliteID2 %in% res_TargHbA1cAdj_Age$Metabolite)))
relat_UntarMetUntarMetAdj <- subset(relat_UntarMetUntarMet, subset=((MetaboliteID1 %in% res_UntargHbA1cAdj_Age$BIOCHEMICAL2)&(MetaboliteID2 %in% res_UntargHbA1cAdj_Age$BIOCHEMICAL2)))
relat_TarMetUntarMetAdj <-  subset(relat_TarMetUntarMet, subset=((MetaboliteID1 %in% res_TargHbA1cAdj_Age$Metabolite)&(MetaboliteID2 %in% res_UntargHbA1cAdj_Age$BIOCHEMICAL2)))

##	MULTIPLE CORRECTION (Qvalue)
relat_ExprsExprsAdj$Exprs_Q <- qvalue(as.numeric(relat_ExprsExprsAdj[,"Exprs_pvalue"]))$qvalues
relat_ProtProtAdj$Prot_Q <- qvalue(as.numeric(relat_ProtProtAdj[,"Prot_pvalue"]))$qvalues
relat_TarMetTarMetAdj$TarMet_Q <- qvalue(as.numeric(relat_TarMetTarMetAdj[,"TarMet_pvalue"]), pi0=1)$qvalues
relat_UntarMetUntarMetAdj$UntarMet_Q <- qvalue(as.numeric(relat_UntarMetUntarMetAdj[,"UntarMet_pvalue"]))$qvalues
relat_TarMetUntarMetAdj$UntarMet_Q <- qvalue(as.numeric(relat_TarMetUntarMetAdj[,"UntarMet_pvalue"]))$qvalues
#relat_ExprsExprsAdj$Exprs_FDR <- p.adjust(as.numeric(relat_ExprsExprsAdj$Exprs_pvalue), method = "BH")
relat_ProtProtAdj$Prot_FDR <- p.adjust(as.numeric(relat_ProtProtAdj$Prot_pvalue), method = "BH")
relat_TarMetTarMetAdj$TarMet_FDR <- p.adjust(as.numeric(relat_TarMetTarMetAdj$TarMet_pvalue), method = "BH")
relat_UntarMetUntarMetAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_UntarMetUntarMetAdj$UntarMet_pvalue), method = "BH")
relat_TarMetUntarMetAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_TarMetUntarMetAdj$UntarMet_pvalue), method = "BH")

sum(relat_ExprsExprsAdj$Exprs_Q < 0.05)
#	108,000,651
sum(relat_ProtProtAdj$Prot_Q < 0.05)
#	139,129
sum(relat_TarMetTarMetAdj$TarMet_Q < 0.05)
#	13,185
sum(relat_UntarMetUntarMetAdj$UntarMet_Q < 0.05)
#	9,108
sum(relat_TarMetUntarMetAdj$UntarMet_Q < 0.05)
#	0
#sum(relat_ExprsExprsAdj$Exprs_FDR < 0.05)
#	
sum(relat_ProtProtAdj$Prot_FDR < 0.05)
#	131,051
sum(relat_TarMetTarMetAdj$TarMet_FDR < 0.05)
#	13,185
sum(relat_UntarMetUntarMetAdj$UntarMet_FDR < 0.05)
#	7,616
sum(relat_TarMetUntarMetAdj$UntarMet_FDR < 0.05)
#	0

write.table(relat_ExprsExprs, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ExprsExprs_Assoc.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_ProtProtAdj, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ProtProt_Assoc.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_TarMetTarMetAdj, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/TarMetTarMet_Assoc.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_UntarMetUntarMetAdj, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/UntarMetUntarMet_Assoc.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_TarMetUntarMetAdj, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/TarMetUntarMet_Assoc.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

##	SUBSET PAIRS SIGNIFICANT FOR ASSOCIATION
relat_ExprsExprsAdj_sub <- subset(relat_ExprsExprsAdj, subset=(Exprs_Q <= 0.05))
#	108,000,651
relat_ProtProtAdj_sub <- subset(relat_ProtProtAdj, subset=(Prot_FDR <= 0.05))
#	131,051
relat_TarMetTarMetAdj_sub <- subset(relat_TarMetTarMetAdj, subset=(TarMet_FDR <= 0.05))
#	13,185
relat_UntarMetUntarMetAdj_sub <- subset(relat_UntarMetUntarMetAdj, subset=(UntarMet_FDR <= 0.05))
#	7,616
relat_TarMetUntarMetAdj_sub <- subset(relat_TarMetUntarMetAdj, subset=(UntarMet_FDR <= 0.05))
#	0

rownames(relat_ExprsExprsAdj_sub) <- 1:nrow(relat_ExprsExprsAdj_sub)
rownames(relat_ProtProtAdj_sub) <- 1:nrow(relat_ProtProtAdj_sub)
rownames(relat_TarMetTarMetAdj_sub) <- 1:nrow(relat_TarMetTarMetAdj_sub)
rownames(relat_UntarMetUntarMetAdj_sub) <- 1:nrow(relat_UntarMetUntarMetAdj_sub)
#rownames(relat_TarMetUntarMetAdj_sub) <- 1:nrow(relat_TarMetUntarMetAdj_sub)

write.table(relat_ExprsExprsAdj_sub, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ExprsExprs_Assoc_Signif.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_ProtProtAdj_sub, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/ProtProt_Assoc_Signif.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_TarMetTarMetAdj_sub, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/TarMetTarMet_Assoc_Signif.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(relat_UntarMetUntarMetAdj_sub, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/UntarMetUntarMet_Assoc_Signif.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
#write.table(relat_TarMetUntarMetAdj_sub, file = "/users/home/dafmic/CrossSectional/RELATIONS/results_TOTALS/TarMetUntarMet_Assoc_Signif.tab", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)



#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#																									BAYES CriT. INFORMATION SCORE
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
split  -l1542867 ExprsExprs_Assoc_Signif.tab ExprsExprs_Assoc_Signif_
# Separated into 70 files - All 1542867 lines apart from last (1542829)
# Run SCRIPT 70 times, manually changing input and output every time
qsub  /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_Exprs_Exprs.sh
qsub  /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_Exprs_Prot.sh
qsub  /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_Exprs_TarMet.sh
qsub  /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_Exprs_UntarMet.sh
qsub  /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_Prot_Prot.sh
qsub  /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_Prot_Prot2.sh	#x2 [21, 75]
qsub  /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_Prot_TarMet.sh
qsub  /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_TarMet_TarMet.sh
qsub  /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_UntarMet_UntarMet.sh
#*****************************************************************************************************************************************************************************************************************************
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_aa_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_1_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ab_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_2_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ac_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_3_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ad_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_4_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ae_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_5_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_af_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_6_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ag_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_7_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ah_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_8_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ai_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_9_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_aj_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_10_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ak_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_11_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_al_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_12_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_am_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_13_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_an_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_14_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ao_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_15_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ap_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_16_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_aq_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_17_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ar_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_18_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_as_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_19_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_at_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_20_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_au_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_21_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_av_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_22_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_aw_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_23_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ax_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_24_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ay_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_25_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_az_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_26_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ba_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_27_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bb_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_28_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bc_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_29_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bd_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_30_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_be_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_31_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bf_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_32_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bg_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_33_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bh_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_34_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bi_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_35_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bj_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_36_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bk_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_37_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bl_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_38_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bm_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_39_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bn_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_40_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bo_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_41_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bp_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_42_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bq_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_43_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_br_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_44_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bs_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_45_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bt_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_46_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bu_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_47_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bv_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_48_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bw_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_49_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bx_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_50_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_by_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_51_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_bz_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_52_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ca_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_53_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cb_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_54_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cc_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_55_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cd_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_56_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ce_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_57_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cf_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_58_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cg_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_59_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ch_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_60_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ci_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_61_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cj_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_62_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ck_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_63_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cl_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_64_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cm_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_65_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cn_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_66_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_co_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_67_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cp_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_68_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cq_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_69_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_cr_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_70_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
#*****************************************************************************************************************************************************************************************************************************
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_{1..70} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Exprs_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
#*****************************************************************************************************************************************************************************************************************************
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Prot_{1..1000} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_Prot_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_TarMet_{1..100} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_TarMet_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_UntarMet_{1..100} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Exprs_UntarMet_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Prot_Prot_{1..100} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Prot_Prot_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Prot_TarMet_{1..10} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_Prot_TarMet_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_TarMet_TarMet_{1..10} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_TarMet_TarMet_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0
echo "cat /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_UntarMet_UntarMet_{1..10} > /users/home/dafmic/CrossSectional/RELATIONS/BIC/BIC_UntarMet_UntarMet_ALL.txt"| qsub -l nodes=1:ppn=1,mem=10gb,walltime=15:00:0