# THIS FILE CONTAINS THE R CODE USED TO IDENTIFY RELATIONSHIPS (ASSOCIATIONS) BETWEEN MOLECULAR PHENOTYPES 
# AND CONTEXT-DEPENDENT EFFECTS (OF AGE, SEX, BMI, HbA1c) MODULATING THESE RELATIONSHIPS
# FOR GENES, PROTEINS AND METABOLITES (TARGETED & UNTARGETED DATASET)
# USING LINEAR REGRESSION AND INTERACTION MODELS
# AND APPLYING MULTIPLE CORRECTION METHODS
# ALSO INCLUDES CODE FOR CONSTRUCTING BARPLOT AND UPSET PLOTS
# FINALLY, THE CODE FOR INTERACTION PLOTS IS INCLUDED

#-------------------------------------------------------------------	LIBRARIES/PACKAGES USED	-------------------------------------------------------------------------
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(UpSetR)
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#	CALCULATE RESIDUALS [remove effects of TECHNICAL covariates only ]

#	1. Gene Expression
resid_genes_Tech <- matrix(0,ncol=3029,nrow=16209)
colnames(resid_genes_Tech) <- c(colnames(anno_genes)[4:5],colnames(raw_genesE))
for(i in 1:nrow(resid_genes_Tech)){
	exprs <- as.matrix(as.double(raw_genesE[i,]))
	lmTemp<- lm(exprs ~ covar_Exprsf$WP + covar_Exprsf$Centre + covar_Exprsf$GC_mean + covar_Exprsf$InsertSizeMode + covar_Exprsf$Rin + covar_Exprsf$Date)
	resid_genes_Tech[i,] <- c(as.character(anno_genes[i,4:5]), as.numeric(resid(lmTemp)))
}
#	2. Proteins
panel <- rbind(c('Olink_CARDIOMETABOLIC','C-Met','CAM'), c('Olink_CARDIOVASCULAR_II','CVD_II','CVD2'),
               c('Olink_CARDIOVASCULAR_III','CVD_III','CVD3'), c('Olink_DEVELOPMENT','Dev','DEV'),
               c('Olink_METABOLISM','Met','MET'))
resid_proteins_Tech <- matrix(0,ncol=3025,nrow=373)
colnames(resid_proteins_Tech) <- c(colnames(anno_proteins)[c(4,6)],colnames(raw_proteins_Imp))
for(i in 1:nrow(resid_proteins_Tech)){
	plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[i,4],4],3], sep='')]))
	prot <- as.matrix(as.double(raw_proteins_Imp[i,]))
	lmTemp<- lm(prot ~ covar_Protf$WP + covar_Protf$Center + plate)
	resid_proteins_Tech[i,] <- c(as.character(anno_proteins[i,c(4,6)]), as.numeric(resid(lmTemp)))
}
#	3. Targeted Metabolites
resid_tarMet_Tech <- matrix(0,ncol=3028,nrow=116)
colnames(resid_tarMet_Tech) <- c(colnames(anno_met_targ),colnames(raw_met_targE))
for(i in 1:nrow(resid_tarMet_Tech)){
	tar <- as.matrix(as.double(raw_met_targE[i,]))
	lmTemp<- lm(tar ~ covar_Targf$WP + covar_Targf$Centre + covar_Targf$Plate.Bar.Code + covar_Targf$date)
	resid_tarMet_Tech[i,] <- c(as.character(anno_met_targ[i,]), as.numeric(resid(lmTemp)))
}
#	4. Untargeted Metabolites
resid_untarMet_Tech <- matrix(0,ncol=2998,nrow=233)
colnames(resid_untarMet_Tech) <- c(colnames(anno_met_untarg)[3:4],colnames(raw_met_untarg_Imp))
for(i in 1:nrow(resid_untarMet_Tech)){
	untar <- as.matrix(as.double(raw_met_untarg_Imp[i,]))
	lmTemp<- lm(untar ~ covar_Untargf$WP + covar_Untargf$Centre + covar_Untargf$Run_Day)
	resid_untarMet_Tech[i,] <- c(as.character(anno_met_untarg[i,3:4]), as.numeric(resid(lmTemp)))
}


#	SUBSET DATASETS ACROSS MOLPHE [INDIVIDUALS WITH BOTH MOLPHE TYPES MEASUREMENTS]

#	GeneExpression:				16,209 x 3,027
#	Protein Levels:					373 x 3,023
#	Targeted Metabolites:		116 x 3,027
#	Untargeted Metabolites:	233 x 2,996

# SUBSET RESIDUALISED EXPRS FOR 3,023 INDIVIDUALS WITH PROTEOMICS AND 2,9996 WITH UNTARGETED METABOLOMICS
resid_genes_Tech_sub3023 <- data.frame(resid_genes_Tech[,1:2],resid_genes_Tech[,c(which(colnames(resid_genes_Tech) %in% colnames(raw_proteins_Imp)))])
colnames(resid_genes_Tech_sub3023)[3:3025] <- colnames(raw_proteins_Imp)
resid_genes_Tech_sub2996 <- data.frame(resid_genes_Tech[,1:2],resid_genes_Tech[,c(which(colnames(resid_genes_Tech) %in% colnames(raw_met_untarg_Imp)))])
colnames(resid_genes_Tech_sub2996)[3:2998] <- colnames(raw_met_untarg_Imp)
# SUBSET RESIDUALISED TARGETED METABOLITES FOR 3,023 INDIVIDUALS WITH PROTEOMICS
resid_tarMet_Tech_sub3023<- data.frame(resid_tarMet_Tech[,1],resid_tarMet_Tech[,c(which(colnames(resid_tarMet_Tech) %in% colnames(raw_proteins_Imp)))])
colnames(resid_tarMet_Tech_sub3023)[2:3024] <- colnames(raw_proteins_Imp)
colnames(resid_tarMet_Tech_sub3023)[1] <- "Metabolite"
# SUBSET RESIDUALISED TARGETED METABOLITES FOR 2,9996 INDIVIDUALS WITH PROTEOMICS
resid_tarMet_Tech_sub2996<- data.frame(resid_tarMet_Tech[,1],resid_tarMet_Tech[,c(which(colnames(resid_tarMet_Tech) %in% colnames(raw_met_untarg_Imp)))])
colnames(resid_tarMet_Tech_sub2996)[2:2997] <- colnames(raw_met_untarg_Imp)
colnames(resid_tarMet_Tech_sub2996)[1] <- "Metabolite"
# SUBSET RESIDUALISED PROTEINS FOR 2,9996  INDIVIDUALS WITH UNTARGETED METABOLOMICS
length(intersect(colnames(resid_untarMet_Tech[,3:2998]), colnames(resid_proteins_Tech[,3:3025])))
# 2992 common IDs
resid_proteins_Tech_sub2992 <- cbind(resid_proteins_Tech[,1:2],resid_proteins_Tech[,c(which(colnames(resid_proteins_Tech) %in% colnames(raw_met_untarg_Imp)))])
colnames(resid_proteins_Tech_sub2992)[3:2994] <- colnames(raw_met_untarg_Imp)
overlap <- intersect(colnames(resid_untarMet_Tech[,3:2998]),colnames(resid_proteins_Tech_sub2992[,3:2994]))
resid_untarMet_Tech2 <- cbind(resid_untarMet_Tech[,1:2], resid_untarMet_Tech[,overlap])
covar_Untargf2 <- subset(covar_Untargf, subset=(rownames(covar_Untargf)%in%overlap))

#		MOLECULAR PHENOTYPE RELATIONSHIPS (ASSOCIATIONS): LINEAR REGRESSION - FUNCTION
lmMolPheAGEx2 <- function(exprs1, exprs2, exprs1_annot, exprs2_annot, exprs2Name, covar){
results <- matrix(0,ncol=3+ncol(exprs1_annot)+ncol(exprs2_annot),nrow=nrow(exprs1)*nrow(exprs2))
colnames(results) <- c(colnames(exprs1_annot),colnames(exprs2_annot),exprs2Name,paste(exprs2Name,"Tvalue",sep="_"),paste(exprs2Name,"pvalue",sep="_"))
	for(i in 1:nrow(exprs1)){
		print(i)
		k= (i-1) * nrow(exprs2)
		exprsA <- as.matrix(as.double(exprs1[i,]))
		for(j in 1:nrow(exprs2)){
			lmx2<- lm(exprsA ~ as.matrix(as.double(exprs2[j,])) + ., data = covar)
			coeff <- summary(lmx2)$coefficients
			results[k+j,] <- c(as.character(exprs1_annot[i,]), as.character(exprs2_annot[j,]), coeff[2,c(1,3:4)])
		}
	}
	return(results)
}

#		MOLECULAR PHENOTYPE RELATIONSHIPS (ASSOCIATIONS): LINEAR REGRESSION - GET RESULTS
##	RUN		EXPS ~ PROT
identical(colnames(resid_genes_Tech_sub3023[,3:3025]),colnames(resid_proteins_Tech[,3:3025]))
relat_ExprsProt <- lmMolPheAGEx2(exprs1 = resid_genes_Tech_sub3023[,3:3025],
																	exprs2 = resid_proteins_Tech[,3:3025],
																	exprs1_annot = resid_genes_Tech_sub3023[,1:2],
																	exprs2_annot = resid_proteins_Tech[,1:2],
																	covar = covar_Protf[1:4],
																	exprs2Name = "Prot")															
relat_ExprsProt <- as.data.frame(relat_ExprsProt)
# RUN		EXPS ~ MET_TAR
identical(colnames(resid_genes_Tech[,3:3029]),colnames(resid_tarMet_Tech[,2:3028]))
relat_ExprsTarMet <- lmMolPheAGEx2(exprs1 = resid_genes_Tech[,3:3029],
																	exprs2 = resid_tarMet_Tech[,2:3028],
																	exprs1_annot = resid_genes_Tech[,1:2],
																	exprs2_annot = as.data.frame(resid_tarMet_Tech[,1]),
																	covar = covar_Targf[1:4],
																	exprs2Name = "TarMet")																
relat_ExprsTarMet <- as.data.frame(relat_ExprsTarMet)
colnames(relat_ExprsTarMet)[3] <- "MetaboliteID"
# RUN		EXPS ~ MET_UNTAR
identical(colnames(resid_genes_Tech_sub2996[,3:2998]),colnames(resid_untarMet_Tech[,3:2998]))
relat_ExprsUntarMet <- lmMolPheAGEx2(exprs1 = resid_genes_Tech_sub2996[,3:2998],
																	exprs2 = resid_untarMet_Tech[,3:2998],
																	exprs1_annot = resid_genes_Tech_sub2996[,1:2],
																	exprs2_annot = resid_untarMet_Tech[,1:2],
																	covar = covar_Untargf[1:4],
																	exprs2Name = "UntarMet")														
relat_ExprsUntarMet <- as.data.frame(relat_ExprsUntarMet)
# RUN		PROT ~ MET_TAR
identical(colnames(resid_proteins_Tech[,3:3025]),colnames(resid_tarMet_Tech_sub3023[,2:3024]))
relat_ProtTarMet <- lmMolPheAGEx2(exprs1 = resid_proteins_Tech[,3:3025],
																	exprs2 = resid_tarMet_Tech_sub3023[,2:3024],
																	exprs1_annot = resid_proteins_Tech[,1:2],
																	exprs2_annot = as.data.frame(resid_tarMet_Tech_sub3023[,1]),
																	covar = covar_Protf[1:4],
																	exprs2Name = "TarMet")														
relat_ProtTarMet <- as.data.frame(relat_ProtTarMet)
colnames(relat_ProtTarMet)[3] <- "MetaboliteID"
# RUN		PROT ~ MET_UNTAR
identical(colnames(resid_proteins_Tech_sub2992[,3:2994]),colnames(resid_untarMet_Tech2[,3:2994]))
relat_ProtUntarMet <- lmMolPheAGEx2(exprs1 = resid_proteins_Tech_sub2992[,3:2994],
																	exprs2 = resid_untarMet_Tech2[,3:2994],
																	exprs1_annot = resid_proteins_Tech_sub2992[,1:2],
																	exprs2_annot = resid_untarMet_Tech2[,1:2],
																	covar = covar_Untargf2[1:4],
																	exprs2Name = "UntarMet")						
relat_ProtUntarMet <- as.data.frame(relat_ProtUntarMet)



#		CONTEXT-DEPENDET EFFECT ON MOLECULAR PHENOTYPE RELATIONSHIPS (ASSOCIATIONS): INTERACTION MODEL - FUNCTION
lmMolPheAGEInter <- function(exprs1, exprs2, exprs1_annot, exprs2_annot, covar, cov_binary = FALSE, biolCovar, exprs2Name, covName){
results <- matrix(0,ncol=9+ncol(exprs1_annot)+ncol(exprs2_annot),nrow=nrow(exprs1)*nrow(exprs2))
colnames(results) <- c(colnames(exprs1_annot),colnames(exprs2_annot),exprs2Name,covName,paste(exprs2Name,covName,sep="*"),paste(exprs2Name,"Tvalue",sep="_"),paste(covName,"Tvalue",sep="_"),paste(paste(exprs2Name,covName,sep="*"),"Tvalue",sep="_"),paste(exprs2Name,"pvalue",sep="_"),paste(covName,"pvalue",sep="_"),paste(paste(exprs2Name,covName,sep="*"),"pvalue",sep="_"))

if(cov_binary == F){
	for(i in 1:nrow(exprs1)){
		print(i)
		k= (i-1) * nrow(exprs2)
		exprsA <- as.matrix(as.double(exprs1[i,]))
		for(j in 1:nrow(exprs2)){
			lmInt <- lm(exprsA ~ scale(as.matrix(as.double(exprs2[j,]))) * scale(covar, scale = FALSE) + ., data = biolCovar)
			coeff <- summary(lmInt)$coefficients
			results[k+j,] <- c(as.character(exprs1_annot[i,]), as.character(exprs2_annot[j,]), coeff[2:4,c(1,3:4)])
		}
	}
}	
else{
	for(i in 1:nrow(exprs1)){
		print(i)
		k= (i-1) * nrow(exprs2)
		exprsA <- as.matrix(as.double(exprs1[i,]))
		for(j in 1:nrow(exprs2)){
			lmInt <- lm(exprsA ~ scale(as.matrix(as.double(exprs2[j,]))) * covar + ., data = biolCovar)
			coeff <- summary(lmInt)$coefficients
			results[k+j,] <- c(as.character(exprs1_annot[i,]), as.character(exprs2_annot[j,]), coeff[2:4,c(1,3:4)])
		}
	}
}
return(results)
}


#		CONTEXT-DEPENDET EFFECT ON MOLECULAR PHENOTYPE RELATIONSHIPS (ASSOCIATIONS): INTERACTION MODEL - GET RESULTS (AGE)
##	RUN		EXPS ~ AGE * PROT
relat_ExprsProtAge <- lmMolPheAGEInter(exprs1 = resid_genes_Tech_sub3023[,3:3025],
																	exprs2 = resid_proteins_Tech[,3:3025],
																	exprs1_annot = resid_genes_Tech_sub3023[,1:2],
																	exprs2_annot = resid_proteins_Tech[,1:2],
																	covar = covar_Protf$Age,
																	biolCovar = covar_Protf[2:4],
																	exprs2Name = "Prot",
																	covName = "Age")																	
relat_ExprsProtAge <- as.data.frame(relat_ExprsProtAge)
# RUN		EXPS ~ AGE * MET_TAR
relat_ExprsTarMetAge <- lmMolPheAGEInter(exprs1 = resid_genes_Tech[,3:3029],
																	exprs2 = resid_tarMet_Tech[,2:3028],
																	exprs1_annot = resid_genes_Tech[,1:2],
																	exprs2_annot = as.data.frame(resid_tarMet_Tech[,1]),
																	covar = covar_Targf$Age,
																	biolCovar = covar_Targf[2:4],
																	exprs2Name = "TarMet",
																	covName = "Age")																	
relat_ExprsTarMetAge <- as.data.frame(relat_ExprsTarMetAge)
colnames(relat_ExprsTarMetAge)[3] <- "MetaboliteID"
# RUN		EXPS ~ AGE * MET_UNTAR
relat_ExprsUntarMetAge <- lmMolPheAGEInter(exprs1 = resid_genes_Tech_sub2996[,3:2998],
																	exprs2 = resid_untarMet_Tech[,3:2998],
																	exprs1_annot = resid_genes_Tech_sub2996[,1:2],
																	exprs2_annot = resid_untarMet_Tech[,1:2],
																	covar = covar_Untargf$Age,
																	biolCovar = covar_Untargf[2:4],
																	exprs2Name = "UntarMet",
																	covName = "Age")																	
relat_ExprsUntarMetAge <- as.data.frame(relat_ExprsUntarMetAge)
# RUN		PROT ~ AGE * MET_TAR
relat_ProtTarMetAge <- lmMolPheAGEInter(exprs1 = resid_proteins_Tech[,3:3025],
																	exprs2 = resid_tarMet_Tech_sub3023[,2:3024],
																	exprs1_annot = resid_proteins_Tech[,1:2],
																	exprs2_annot = as.data.frame(resid_tarMet_Tech_sub3023[,1]),
																	covar = covar_Protf$Age,
																	biolCovar = covar_Protf[2:4],
																	exprs2Name = "TarMet",
																	covName = "Age")																	
relat_ProtTarMetAge <- as.data.frame(relat_ProtTarMetAge)
colnames(relat_ProtTarMetAge)[3] <- "MetaboliteID"
# RUN		PROT ~ AGE * MET_UNTAR
relat_ProtUntarMetAge <- lmMolPheAGEInter(exprs1 = resid_proteins_Tech_sub2992[,3:2994],
																	exprs2 = resid_untarMet_Tech2[,3:2994],
																	exprs1_annot = resid_proteins_Tech_sub2992[,1:2],
																	exprs2_annot = resid_untarMet_Tech2[,1:2],
																	covar = covar_Untargf2$Age,
																	biolCovar = covar_Untargf2[2:4],
																	exprs2Name = "UntarMet",
																	covName = "Age")																	
relat_ProtUntarMetAge <- as.data.frame(relat_ProtUntarMetAge)


#		CONTEXT-DEPENDET EFFECT ON MOLECULAR PHENOTYPE RELATIONSHIPS (ASSOCIATIONS): INTERACTION MODEL - GET RESULTS (SEX)
##	RUN		EXPS ~ SEX * PROT
relat_ExprsProtSex <- lmMolPheAGEInter(exprs1 = resid_genes_Tech_sub3023[,3:3025],
																	exprs2 = resid_proteins_Tech[,3:3025],
																	exprs1_annot = resid_genes_Tech_sub3023[,1:2],
																	exprs2_annot = resid_proteins_Tech[,1:2],
																	covar = covar_Protf$Sex,
																	cov_binary = TRUE,
																	biolCovar = covar_Protf[c(1,3:4)],
																	exprs2Name = "Prot",
																	covName = "Sex")																	
relat_ExprsProtSex <- as.data.frame(relat_ExprsProtSex)
# RUN		EXPS ~ SEX * MET_TAR
relat_ExprsTarMetSex <- lmMolPheAGEInter(exprs1 = resid_genes_Tech[,3:3029],
																	exprs2 = resid_tarMet_Tech[,2:3028],
																	exprs1_annot = resid_genes_Tech[,1:2],
																	exprs2_annot = as.data.frame(resid_tarMet_Tech[,1]),
																	covar = covar_Targf$Sex,
																	cov_binary = TRUE,
																	biolCovar = covar_Targf[c(1,3:4)],
																	exprs2Name = "TarMet",
																	covName = "Sex")																	
relat_ExprsTarMetSex <- as.data.frame(relat_ExprsTarMetSex)
colnames(relat_ExprsTarMetSex)[3] <- "MetaboliteID"
# RUN		EXPS ~ SEX * MET_UNTAR
relat_ExprsUntarMetSex <- lmMolPheAGEInter(exprs1 = resid_genes_Tech_sub2996[,3:2998],
																	exprs2 = resid_untarMet_Tech[,3:2998],
																	exprs1_annot = resid_genes_Tech_sub2996[,1:2],
																	exprs2_annot = resid_untarMet_Tech[,1:2],
																	covar = covar_Untargf$Sex,
																	cov_binary = TRUE,
																	biolCovar = covar_Untargf[c(1,3:4)],
																	exprs2Name = "UntarMet",
																	covName = "Sex")																	
relat_ExprsUntarMetSex <- as.data.frame(relat_ExprsUntarMetSex)
# RUN		PROT ~ SEX * MET_TAR
relat_ProtTarMetSex <- lmMolPheAGEInter(exprs1 = resid_proteins_Tech[,3:3025],
																	exprs2 = resid_tarMet_Tech_sub3023[,2:3024],
																	exprs1_annot = resid_proteins_Tech[,1:2],
																	exprs2_annot = as.data.frame(resid_tarMet_Tech_sub3023[,1]),
																	covar = covar_Protf$Sex,
																	cov_binary = TRUE,
																	biolCovar = covar_Protf[c(1,3:4)],
																	exprs2Name = "TarMet",
																	covName = "Sex")																	
relat_ProtTarMetSex <- as.data.frame(relat_ProtTarMetSex)
colnames(relat_ProtTarMetSex)[3] <- "MetaboliteID"
# RUN		PROT ~ SEX * MET_UNTAR
relat_ProtUntarMetSex <- lmMolPheAGEInter(exprs1 = resid_proteins_Tech_sub2992[,3:2994],
																	exprs2 = resid_untarMet_Tech2[,3:2994],
																	exprs1_annot = resid_proteins_Tech_sub2992[,1:2],
																	exprs2_annot = resid_untarMet_Tech2[,1:2],
																	covar = as.factor(covar_Untargf2$Sex),
																	cov_binary = TRUE,
																	biolCovar = covar_Untargf2[c(1,3:4)],
																	exprs2Name = "UntarMet",
																	covName = "Sex")																	
relat_ProtUntarMetSex <- as.data.frame(relat_ProtUntarMetSex)



#		CONTEXT-DEPENDET EFFECT ON MOLECULAR PHENOTYPE RELATIONSHIPS (ASSOCIATIONS): INTERACTION MODEL - GET RESULTS (BMI)
##	RUN		EXPS ~ BMI * PROT
relat_ExprsProtBMI <- lmMolPheAGEInter(exprs1 = resid_genes_Tech_sub3023[,3:3025],
																	exprs2 = resid_proteins_Tech[,3:3025],
																	exprs1_annot = resid_genes_Tech_sub3023[,1:2],
																	exprs2_annot = resid_proteins_Tech[,1:2],
																	covar = covar_Protf$BMI,
																	biolCovar = covar_Protf[c(1:2,4)],
																	exprs2Name = "Prot",
																	covName = "BMI")																	
relat_ExprsProtBMI <- as.data.frame(relat_ExprsProtBMI)
# RUN		EXPS ~ BMI * MET_TAR
relat_ExprsTarMetBMI <- lmMolPheAGEInter(exprs1 = resid_genes_Tech[,3:3029],
																	exprs2 = resid_tarMet_Tech[,2:3028],
																	exprs1_annot = resid_genes_Tech[,1:2],
																	exprs2_annot = as.data.frame(resid_tarMet_Tech[,1]),
																	covar = covar_Targf$BMI,
																	biolCovar = covar_Targf[c(1:2,4)],
																	exprs2Name = "TarMet",
																	covName = "BMI")																	
relat_ExprsTarMetBMI <- as.data.frame(relat_ExprsTarMetBMI)
colnames(relat_ExprsTarMetBMI)[3] <- "MetaboliteID"
# RUN		EXPS ~ BMI * MET_UNTAR
relat_ExprsUntarMetBMI <- lmMolPheAGEInter(exprs1 = resid_genes_Tech_sub2996[,3:2998],
																	exprs2 = resid_untarMet_Tech[,3:2998],
																	exprs1_annot = resid_genes_Tech_sub2996[,1:2],
																	exprs2_annot = resid_untarMet_Tech[,1:2],
																	covar = covar_Untargf$BMI,
																	biolCovar = covar_Untargf[c(1:2,4)],
																	exprs2Name = "UntarMet",
																	covName = "BMI")																	
relat_ExprsUntarMetBMI <- as.data.frame(relat_ExprsUntarMetBMI)
# RUN		PROT ~ BMI * MET_TAR
relat_ProtTarMetBMI <- lmMolPheAGEInter(exprs1 = resid_proteins_Tech[,3:3025],
																	exprs2 = resid_tarMet_Tech_sub3023[,2:3024],
																	exprs1_annot = resid_proteins_Tech[,1:2],
																	exprs2_annot = as.data.frame(resid_tarMet_Tech_sub3023[,1]),
																	covar = covar_Protf$BMI,
																	biolCovar = covar_Protf[c(1:2,4)],
																	exprs2Name = "TarMet",
																	covName = "BMI")																	
relat_ProtTarMetBMI <- as.data.frame(relat_ProtTarMetBMI)
colnames(relat_ProtTarMetBMI)[3] <- "MetaboliteID"
# RUN		PROT ~ BMI * MET_UNTAR
relat_ProtUntarMetBMI <- lmMolPheAGEInter(exprs1 = resid_proteins_Tech_sub2992[,3:2994],
																	exprs2 = resid_untarMet_Tech2[,3:2994],
																	exprs1_annot = resid_proteins_Tech_sub2992[,1:2],
																	exprs2_annot = resid_untarMet_Tech2[,1:2],
																	covar = covar_Untargf2$BMI,
																	biolCovar = covar_Untargf2[c(1:2,4)],
																	exprs2Name = "UntarMet",
																	covName = "BMI")																	
relat_ProtUntarMetBMI <- as.data.frame(relat_ProtUntarMetBMI)


#		CONTEXT-DEPENDET EFFECT ON MOLECULAR PHENOTYPE RELATIONSHIPS (ASSOCIATIONS): INTERACTION MODEL - GET RESULTS (HbA1c)
##	RUN		EXPS ~ HbA1c * PROT
relat_ExprsProtHbA1c <- lmMolPheAGEInter(exprs1 = resid_genes_Tech_sub3023[,3:3025],
																	exprs2 = resid_proteins_Tech[,3:3025],
																	exprs1_annot = resid_genes_Tech_sub3023[,1:2],
																	exprs2_annot = resid_proteins_Tech[,1:2],
																	covar = covar_Protf$HbA1c,
																	biolCovar = covar_Protf[1:3],
																	exprs2Name = "Prot",
																	covName = "HbA1c")																	
relat_ExprsProtHbA1c <- as.data.frame(relat_ExprsProtHbA1c)
# RUN		EXPS ~ HbA1c * MET_TAR
relat_ExprsTarMetHbA1c <- lmMolPheAGEInter(exprs1 = resid_genes_Tech[,3:3029],
																	exprs2 = resid_tarMet_Tech[,2:3028],
																	exprs1_annot = resid_genes_Tech[,1:2],
																	exprs2_annot = as.data.frame(resid_tarMet_Tech[,1]),
																	covar = covar_Targf$HbA1c,
																	biolCovar = covar_Targf[1:3],
																	exprs2Name = "TarMet",
																	covName = "HbA1c")																	
relat_ExprsTarMetHbA1c <- as.data.frame(relat_ExprsTarMetHbA1c)
colnames(relat_ExprsTarMetHbA1c)[3] <- "MetaboliteID"
# RUN		EXPS ~ HbA1c * MET_UNTAR
relat_ExprsUntarMetHbA1c <- lmMolPheAGEInter(exprs1 = resid_genes_Tech_sub2996[,3:2998],
																	exprs2 = resid_untarMet_Tech[,3:2998],
																	exprs1_annot = resid_genes_Tech_sub2996[,1:2],
																	exprs2_annot = resid_untarMet_Tech[,1:2],
																	covar = covar_Untargf$HbA1c,
																	biolCovar = covar_Untargf[1:3],
																	exprs2Name = "UntarMet",
																	covName = "HbA1c")																	
relat_ExprsUntarMetHbA1c <- as.data.frame(relat_ExprsUntarMetHbA1c)
# RUN		PROT ~ HbA1c * MET_TAR
relat_ProtTarMetHbA1c <- lmMolPheAGEInter(exprs1 = resid_proteins_Tech[,3:3025],
																	exprs2 = resid_tarMet_Tech_sub3023[,2:3024],
																	exprs1_annot = resid_proteins_Tech[,1:2],
																	exprs2_annot = as.data.frame(resid_tarMet_Tech_sub3023[,1]),
																	covar = covar_Protf$HbA1c,
																	biolCovar = covar_Protf[1:3],
																	exprs2Name = "TarMet",
																	covName = "HbA1c")																	
relat_ProtTarMetHbA1c <- as.data.frame(relat_ProtTarMetHbA1c)
colnames(relat_ProtTarMetHbA1c)[3] <- "MetaboliteID"
# RUN		PROT ~ HbA1c * MET_UNTAR
relat_ProtUntarMetHbA1c <- lmMolPheAGEInter(exprs1 = resid_proteins_Tech_sub2992[,3:2994],
																	exprs2 = resid_untarMet_Tech2[,3:2994],
																	exprs1_annot = resid_proteins_Tech_sub2992[,1:2],
																	exprs2_annot = resid_untarMet_Tech2[,1:2],
																	covar = covar_Untargf2$HbA1c,
																	biolCovar = covar_Untargf2[1:3],
																	exprs2Name = "UntarMet",
																	covName = "HbA1c")																	
relat_ProtUntarMetHbA1c <- as.data.frame(relat_ProtUntarMetHbA1c)



#	MULTIPLE CORRECTION
relat_ExprsProt$Prot_FDR <- p.adjust(as.numeric(relat_ExprsProt$Prot_pvalue), method = "BH")
relat_ExprsTarMet$TarMet_FDR <- p.adjust(as.numeric(relat_ExprsTarMet$TarMet_pvalue), method = "BH")
relat_ExprsUntarMet$UntarMet_FDR <- p.adjust(as.numeric(relat_ExprsUntarMet$UntarMet_pvalue), method = "BH")
relat_ProtTarMet$TarMet_FDR <- p.adjust(as.numeric(relat_ProtTarMet$TarMet_pvalue), method = "BH")
relat_ProtUntarMet$UntarMet_FDR <- p.adjust(as.numeric(relat_ProtUntarMet$UntarMet_pvalue), method = "BH")
relat_ExprsProt$Prot_Q <- qvalue(as.numeric(unlist(relat_ExprsProt[,"Prot_pvalue"])))$qvalues
relat_ExprsTarMet$TarMet_Q <- qvalue(as.numeric(unlist(relat_ExprsTarMet[,"TarMet_pvalue"])))$qvalues
relat_ExprsUntarMet$UntarMet_Q <- qvalue(as.numeric(unlist(relat_ExprsUntarMet[,"UntarMet_pvalue"])))$qvalues
relat_ProtTarMet$TarMet_Q <- qvalue(as.numeric(unlist(relat_ProtTarMet[,"TarMet_pvalue"])))$qvalues
relat_ProtUntarMet$UntarMet_Q <- qvalue(as.numeric(unlist(relat_ProtUntarMet[,"UntarMet_pvalue"])))$qvalues
#relat_ExprsExprs$Exprs_FDR <- p.adjust(as.numeric(relat_ExprsExprs$Exprs_pvalue), method = "BH")
relat_ProtProt$Prot_FDR <- p.adjust(as.numeric(relat_ProtProt$Prot_pvalue), method = "BH")
relat_TarMetTarMet$TarMet_FDR <- p.adjust(as.numeric(relat_TarMetTarMet$TarMet_pvalue), method = "BH")
relat_UntarMetUntarMet$UntarMet_FDR <- p.adjust(as.numeric(relat_UntarMetUntarMet$UntarMet_pvalue), method = "BH")
relat_TarMetUntarMet$UntarMet_FDR <- p.adjust(as.numeric(relat_TarMetUntarMet$UntarMet_pvalue), method = "BH")
#relat_ExprsExprs$Exprs_Q <- qvalue(as.numeric(unlist(relat_ExprsExprs[,"Exprs_pvalue"])))$qvalues
relat_ProtProt$Prot_Q <- qvalue(as.numeric(unlist(relat_ProtProt[,"Prot_pvalue"])))$qvalues
relat_TarMetTarMet$TarMet_Q <- qvalue(as.numeric(unlist(relat_TarMetTarMet[,"TarMet_pvalue"])), pi0=1)$qvalues
relat_UntarMetUntarMet$UntarMet_Q <- qvalue(as.numeric(unlist(relat_UntarMetUntarMet[,"UntarMet_pvalue"])))$qvalues
relat_TarMetUntarMet$UntarMet_Q <- qvalue(as.numeric(unlist(relat_TarMetUntarMet[,"UntarMet_pvalue"])))$qvalues

sum(relat_ExprsProt$Prot_Q < 0.05)
#	1,346,330/6,045,957					22.27%
sum(relat_ExprsTarMet$TarMet_Q < 0.05)
#	293,412/1,880,244	
sum(relat_ExprsUntarMet$UntarMet_Q < 0.05)
#	276,705/3,776,697	
# Exprs ~ Metab (merged): 570,117/5,656,941	10.08%
sum(relat_ProtTarMet$TarMet_Q < 0.05)
#	25,349/43,268
sum(relat_ProtUntarMet$UntarMet_Q < 0.05)
#	2/86,909
#Prot ~ Metab (merged): 25,351/130,177	19.47%
#------------------------------------------------------
sum(relat_ExprsProt$Prot_FDR < 0.05)
#	1,098,178/6,045,957       18.16%
sum(relat_ExprsTarMet$TarMet_FDR < 0.05)
#	227,634/1,880,244         12.11%
sum(relat_ExprsUntarMet$UntarMet_FDR < 0.05)
#	227,909/3,776,697         6.03%
# Exprs ~ Metab (merged): 455,543/5,656,941	8.05%
sum(relat_ProtTarMet$TarMet_FDR < 0.05)
#	20,056/43,268             46.35%
sum(relat_ProtUntarMet$UntarMet_FDR < 0.05)
#	2/86,909                  0%
#Prot ~ Metab (merged): 20,058/130,177	15.41%

MOLECULES2 <- c(as.vector(unlist(subset(relat_ExprsProt, select=(GeneID), subset=(Prot_FDR < 0.05))[,1])),
               as.vector(unlist(subset(relat_ExprsProt, select=(ProteinID), subset=(Prot_FDR < 0.05))[,1])),
               as.vector(unlist(subset(relat_ExprsTarMet, select=(GeneID), subset=(TarMet_FDR < 0.05))[,1])),
               as.vector(unlist(subset(relat_ExprsTarMet, select=(MetaboliteID), subset=(TarMet_FDR < 0.05))[,1])),
               as.vector(unlist(subset(relat_ExprsUntarMet, select=(GeneID), subset=(UntarMet_FDR < 0.05))[,1])),
               as.vector(unlist(subset(relat_ExprsUntarMet, select=(BIOCHEMICAL2), subset=(UntarMet_FDR < 0.05))[,1])),
               as.vector(unlist(subset(relat_ProtTarMet, select=(MetaboliteID), subset=(TarMet_FDR < 0.05))[,1])),
               as.vector(unlist(subset(relat_ProtTarMet, select=(ProteinID), subset=(TarMet_FDR < 0.05))[,1])),
               as.vector(unlist(subset(relat_ProtUntarMet, select=(BIOCHEMICAL2), subset=(UntarMet_FDR < 0.05))[,1])),
               as.vector(unlist(subset(relat_ProtUntarMet, select=(ProteinID), subset=(UntarMet_FDR < 0.05))[,1])))
MOLECULES2 <- unique(MOLECULES2)# 16,931 unique molecules

# TOTAL MOLECULES: 16,931
# UNIQUE MOLECULES IN SIGNIF PAIRS: 16,931
# 100 %
# UNIQUE MOLECULES IN SIGNIF PAIRS (FDR): 16,927
# 100 %

GP <- subset(relat_ExprsProt, subset=(Prot_Q < 0.05))
GP <- subset(GP, subset=(GeneName==ProteinName))
# 153 Gene-Protein PAIRS from ExprProt

GP2 <- subset(relat_ExprsProt, subset=(Prot_FDR < 0.05))
GP2 <- subset(GP2, subset=(GeneName==ProteinName))
# 139 Gene-Protein PAIRS from ExprProt

relat_ExprsProtAgeAdj <- relat_ExprsProtAge
relat_ExprsProtAgeAdj$Prot_Q <- qvalue(as.numeric(unlist(relat_ExprsProtAgeAdj[,"Prot_pvalue"])))$qvalues
relat_ExprsProtAgeAdj$Age_Q <- qvalue(as.numeric(unlist(relat_ExprsProtAgeAdj[,"Age_pvalue"])))$qvalues
relat_ExprsProtAgeAdj$"Prot*Age_Q" <- qvalue(as.numeric(unlist(relat_ExprsProtAgeAdj[,"Prot*Age_pvalue"])))$qvalues
relat_ExprsProtAgeAdj$Prot_FDR <- p.adjust(as.numeric(relat_ExprsProtAgeAdj$Prot_pvalue), method = "BH")
relat_ExprsProtAgeAdj$Age_FDR <- p.adjust(as.numeric(relat_ExprsProtAgeAdj$Age_pvalue), method = "BH")
relat_ExprsProtAgeAdj$"Prot*Age_FDR" <- p.adjust(as.numeric(relat_ExprsProtAgeAdj$"Prot*Age_pvalue"), method = "BH")

relat_ExprsTarMetAgeAdj <- relat_ExprsTarMetAge
relat_ExprsTarMetAgeAdj$TarMet_Q <- qvalue(as.numeric(unlist(relat_ExprsTarMetAgeAdj[,"TarMet_pvalue"])))$qvalues
relat_ExprsTarMetAgeAdj$Age_Q <- qvalue(as.numeric(unlist(relat_ExprsTarMetAgeAdj[,"Age_pvalue"])))$qvalues
relat_ExprsTarMetAgeAdj$"TarMet*Age_Q" <- qvalue(as.numeric(unlist(relat_ExprsTarMetAgeAdj[,"TarMet*Age_pvalue"])))$qvalues
relat_ExprsTarMetAgeAdj$TarMet_FDR <- p.adjust(as.numeric(relat_ExprsTarMetAgeAdj$"TarMet_pvalue"), method = "BH")
relat_ExprsTarMetAgeAdj$Age_FDR <- p.adjust(as.numeric(relat_ExprsTarMetAgeAdj$"Age_pvalue"), method = "BH")
relat_ExprsTarMetAgeAdj$"TarMet*Age_FDR" <- p.adjust(as.numeric(relat_ExprsTarMetAgeAdj$"TarMet*Age_pvalue"), method = "BH")

relat_ExprsUntarMetAgeAdj <- relat_ExprsUntarMetAge
relat_ExprsUntarMetAgeAdj$UntarMet_Q <- qvalue(as.numeric(unlist(relat_ExprsUntarMetAgeAdj[,"UntarMet_pvalue"])))$qvalues
relat_ExprsUntarMetAgeAdj$Age_Q <- qvalue(as.numeric(unlist(relat_ExprsUntarMetAgeAdj[,"Age_pvalue"])))$qvalues
relat_ExprsUntarMetAgeAdj$"UntarMet*Age_Q" <- qvalue(as.numeric(unlist(relat_ExprsUntarMetAgeAdj[,"UntarMet*Age_pvalue"])))$qvalues
relat_ExprsUntarMetAgeAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_ExprsUntarMetAgeAdj$"UntarMet_pvalue"), method = "BH")
relat_ExprsUntarMetAgeAdj$Age_FDR <- p.adjust(as.numeric(relat_ExprsUntarMetAgeAdj$"Age_pvalue"), method = "BH")
relat_ExprsUntarMetAgeAdj$"UntarMet*Age_FDR" <- p.adjust(as.numeric(relat_ExprsUntarMetAgeAdj$"UntarMet*Age_pvalue"), method = "BH")

relat_ProtTarMetAgeAdj <- relat_ProtTarMetAge
relat_ProtTarMetAgeAdj$TarMet_Q <- qvalue(as.numeric(unlist(relat_ProtTarMetAgeAdj[,"TarMet_pvalue"])))$qvalues
relat_ProtTarMetAgeAdj$Age_Q <- qvalue(as.numeric(unlist(relat_ProtTarMetAgeAdj[,"Age_pvalue"])))$qvalues
relat_ProtTarMetAgeAdj$"TarMet*Age_Q" <- qvalue(as.numeric(unlist(relat_ProtTarMetAgeAdj[,"TarMet*Age_pvalue"])))$qvalues
relat_ProtTarMetAgeAdj$TarMet_FDR <- p.adjust(as.numeric(relat_ProtTarMetAgeAdj$"TarMet_pvalue"), method = "BH")
relat_ProtTarMetAgeAdj$Age_FDR <- p.adjust(as.numeric(relat_ProtTarMetAgeAdj$"Age_pvalue"), method = "BH")
relat_ProtTarMetAgeAdj$"TarMet*Age_FDR" <- p.adjust(as.numeric(relat_ProtTarMetAgeAdj$"TarMet*Age_pvalue"), method = "BH")

relat_ProtUntarMetAgeAdj <- relat_ProtUntarMetAge
relat_ProtUntarMetAgeAdj$UntarMet_Q <- qvalue(as.numeric(unlist(relat_ProtUntarMetAgeAdj[,"UntarMet_pvalue"])))$qvalues
relat_ProtUntarMetAgeAdj$Age_Q <- qvalue(as.numeric(unlist(relat_ProtUntarMetAgeAdj[,"Age_pvalue"])))$qvalues
relat_ProtUntarMetAgeAdj$"UntarMet*Age_Q" <- qvalue(as.numeric(unlist(relat_ProtUntarMetAgeAdj[,"UntarMet*Age_pvalue"])))$qvalues
relat_ProtUntarMetAgeAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_ProtUntarMetAgeAdj$"UntarMet_pvalue"), method = "BH")
relat_ProtUntarMetAgeAdj$Age_FDR <- p.adjust(as.numeric(relat_ProtUntarMetAgeAdj$"Age_pvalue"), method = "BH")
relat_ProtUntarMetAgeAdj$"UntarMet*Age_FDR" <- p.adjust(as.numeric(relat_ProtUntarMetAgeAdj$"UntarMet*Age_pvalue"), method = "BH")

relat_ProtProtAgeAdj <- relat_ProtProtAge
relat_ProtProtAgeAdj$Prot_Q <- qvalue(as.numeric(unlist(relat_ProtProtAgeAdj[,"Prot_pvalue"])))$qvalues
relat_ProtProtAgeAdj$Age_Q <- qvalue(as.numeric(unlist(relat_ProtProtAgeAdj[,"Age_pvalue"])))$qvalues
relat_ProtProtAgeAdj$"Prot*Age_Q" <- qvalue(as.numeric(unlist(relat_ProtProtAgeAdj[,"Prot*Age_pvalue"])))$qvalues
relat_ProtProtAgeAdj$Prot_FDR <- p.adjust(as.numeric(relat_ProtProtAgeAdj$Prot_pvalue), method = "BH")
relat_ProtProtAgeAdj$Age_FDR <- p.adjust(as.numeric(relat_ProtProtAgeAdj$Age_pvalue), method = "BH")
relat_ProtProtAgeAdj$"Prot*Age_FDR" <- p.adjust(as.numeric(relat_ProtProtAgeAdj$"Prot*Age_pvalue"), method = "BH")
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
relat_ExprsProtSexAdj <- relat_ExprsProtSex
relat_ExprsProtSexAdj$Prot_Q <- qvalue(as.numeric(unlist(relat_ExprsProtSexAdj[,"Prot_pvalue"])))$qvalues
relat_ExprsProtSexAdj$Sex_Q <- qvalue(as.numeric(unlist(relat_ExprsProtSexAdj[,"Sex_pvalue"])))$qvalues
relat_ExprsProtSexAdj$"Prot*Sex_Q" <- qvalue(as.numeric(unlist(relat_ExprsProtSexAdj[,"Prot*Sex_pvalue"])))$qvalues
relat_ExprsProtSexAdj$Prot_FDR <- p.adjust(as.numeric(relat_ExprsProtSexAdj$Prot_pvalue), method = "BH")
relat_ExprsProtSexAdj$Sex_FDR <- p.adjust(as.numeric(relat_ExprsProtSexAdj$Sex_pvalue), method = "BH")
relat_ExprsProtSexAdj$"Prot*Sex_FDR" <- p.adjust(as.numeric(relat_ExprsProtSexAdj$"Prot*Sex_pvalue"), method = "BH")

relat_ExprsTarMetSexAdj <- relat_ExprsTarMetSex
relat_ExprsTarMetSexAdj$TarMet_Q <- qvalue(as.numeric(unlist(relat_ExprsTarMetSexAdj[,"TarMet_pvalue"])))$qvalues
relat_ExprsTarMetSexAdj$Sex_Q <- qvalue(as.numeric(unlist(relat_ExprsTarMetSexAdj[,"Sex_pvalue"])))$qvalues
relat_ExprsTarMetSexAdj$"TarMet*Sex_Q" <- qvalue(as.numeric(unlist(relat_ExprsTarMetSexAdj[,"TarMet*Sex_pvalue"])))$qvalues
relat_ExprsTarMetSexAdj$TarMet_FDR <- p.adjust(as.numeric(relat_ExprsTarMetSexAdj$TarMet_pvalue), method = "BH")
relat_ExprsTarMetSexAdj$Sex_FDR <- p.adjust(as.numeric(relat_ExprsTarMetSexAdj$Sex_pvalue), method = "BH")
relat_ExprsTarMetSexAdj$"TarMet*Sex_FDR" <- p.adjust(as.numeric(relat_ExprsTarMetSexAdj$"TarMet*Sex_pvalue"), method = "BH")

relat_ExprsUntarMetSexAdj <- relat_ExprsUntarMetSex
relat_ExprsUntarMetSexAdj$UntarMet_Q <- qvalue(as.numeric(unlist(relat_ExprsUntarMetSexAdj[,"UntarMet_pvalue"])))$qvalues
relat_ExprsUntarMetSexAdj$Sex_Q <- qvalue(as.numeric(unlist(relat_ExprsUntarMetSexAdj[,"Sex_pvalue"])))$qvalues
relat_ExprsUntarMetSexAdj$"UntarMet*Sex_Q" <- qvalue(as.numeric(unlist(relat_ExprsUntarMetSexAdj[,"UntarMet*Sex_pvalue"])))$qvalues
relat_ExprsUntarMetSexAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_ExprsUntarMetSexAdj$UntarMet_pvalue), method = "BH")
relat_ExprsUntarMetSexAdj$Sex_FDR <- p.adjust(as.numeric(relat_ExprsUntarMetSexAdj$Sex_pvalue), method = "BH")
relat_ExprsUntarMetSexAdj$"UntarMet*Sex_FDR" <- p.adjust(as.numeric(relat_ExprsUntarMetSexAdj$"UntarMet*Sex_pvalue"), method = "BH")

relat_ProtTarMetSexAdj <- relat_ProtTarMetSex
relat_ProtTarMetSexAdj$TarMet_Q <- qvalue(as.numeric(unlist(relat_ProtTarMetSexAdj[,"TarMet_pvalue"])))$qvalues
relat_ProtTarMetSexAdj$Sex_Q <- qvalue(as.numeric(unlist(relat_ProtTarMetSexAdj[,"Sex_pvalue"])))$qvalues
relat_ProtTarMetSexAdj$"TarMet*Sex_Q" <- qvalue(as.numeric(unlist(relat_ProtTarMetSexAdj[,"TarMet*Sex_pvalue"])))$qvalues
relat_ProtTarMetSexAdj$TarMet_FDR <- p.adjust(as.numeric(relat_ProtTarMetSexAdj$TarMet_pvalue), method = "BH")
relat_ProtTarMetSexAdj$Sex_FDR <- p.adjust(as.numeric(relat_ProtTarMetSexAdj$Sex_pvalue), method = "BH")
relat_ProtTarMetSexAdj$"TarMet*Sex_FDR" <- p.adjust(as.numeric(relat_ProtTarMetSexAdj$"TarMet*Sex_pvalue"), method = "BH")

relat_ProtUntarMetSexAdj <- relat_ProtUntarMetSex
relat_ProtUntarMetSexAdj$UntarMet_Q <- qvalue(as.numeric(unlist(relat_ProtUntarMetSexAdj[,"UntarMet_pvalue"])))$qvalues
relat_ProtUntarMetSexAdj$Sex_Q <- qvalue(as.numeric(unlist(relat_ProtUntarMetSexAdj[,"Sex_pvalue"])))$qvalues
relat_ProtUntarMetSexAdj$"UntarMet*Sex_Q" <- qvalue(as.numeric(unlist(relat_ProtUntarMetSexAdj[,"UntarMet*Sex_pvalue"])))$qvalues
relat_ProtUntarMetSexAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_ProtUntarMetSexAdj$UntarMet_pvalue), method = "BH")
relat_ProtUntarMetSexAdj$Sex_FDR <- p.adjust(as.numeric(relat_ProtUntarMetSexAdj$Sex_pvalue), method = "BH")
relat_ProtUntarMetSexAdj$"UntarMet*Sex_FDR" <- p.adjust(as.numeric(relat_ProtUntarMetSexAdj$"UntarMet*Sex_pvalue"), method = "BH")
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
relat_ExprsProtBMIAdj <- relat_ExprsProtBMI
relat_ExprsProtBMIAdj$Prot_Q <- qvalue(as.numeric(unlist(relat_ExprsProtBMIAdj[,"Prot_pvalue"])))$qvalues
relat_ExprsProtBMIAdj$BMI_Q <- qvalue(as.numeric(unlist(relat_ExprsProtBMIAdj[,"BMI_pvalue"])))$qvalues
relat_ExprsProtBMIAdj$"Prot*BMI_Q" <- qvalue(as.numeric(unlist(relat_ExprsProtBMIAdj[,"Prot*BMI_pvalue"])))$qvalues
relat_ExprsProtBMIAdj$Prot_FDR <- p.adjust(as.numeric(relat_ExprsProtBMIAdj$Prot_pvalue), method = "BH")
relat_ExprsProtBMIAdj$BMI_FDR <- p.adjust(as.numeric(relat_ExprsProtBMIAdj$BMI_pvalue), method = "BH")
relat_ExprsProtBMIAdj$"Prot*BMI_FDR" <- p.adjust(as.numeric(relat_ExprsProtBMIAdj$"Prot*BMI_pvalue"), method = "BH")

relat_ExprsTarMetBMIAdj <- relat_ExprsTarMetBMI
relat_ExprsTarMetBMIAdj$TarMet_Q <- qvalue(as.numeric(unlist(relat_ExprsTarMetBMIAdj[,"TarMet_pvalue"])))$qvalues
relat_ExprsTarMetBMIAdj$BMI_Q <- qvalue(as.numeric(unlist(relat_ExprsTarMetBMIAdj[,"BMI_pvalue"])))$qvalues
relat_ExprsTarMetBMIAdj$"TarMet*BMI_Q" <- qvalue(as.numeric(unlist(relat_ExprsTarMetBMIAdj[,"TarMet*BMI_pvalue"])))$qvalues
relat_ExprsTarMetBMIAdj$TarMet_FDR <- p.adjust(as.numeric(relat_ExprsTarMetBMIAdj$TarMet_pvalue), method = "BH")
relat_ExprsTarMetBMIAdj$BMI_FDR <- p.adjust(as.numeric(relat_ExprsTarMetBMIAdj$BMI_pvalue), method = "BH")
relat_ExprsTarMetBMIAdj$"TarMet*BMI_FDR" <- p.adjust(as.numeric(relat_ExprsTarMetBMIAdj$"TarMet*BMI_pvalue"), method = "BH")

relat_ExprsUntarMetBMIAdj <- relat_ExprsUntarMetBMI
relat_ExprsUntarMetBMIAdj$UntarMet_Q <- qvalue(as.numeric(unlist(relat_ExprsUntarMetBMIAdj[,"UntarMet_pvalue"])))$qvalues
relat_ExprsUntarMetBMIAdj$BMI_Q <- qvalue(as.numeric(unlist(relat_ExprsUntarMetBMIAdj[,"BMI_pvalue"])))$qvalues
relat_ExprsUntarMetBMIAdj$"UntarMet*BMI_Q" <- qvalue(as.numeric(unlist(relat_ExprsUntarMetBMIAdj[,"UntarMet*BMI_pvalue"])))$qvalues
relat_ExprsUntarMetBMIAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_ExprsUntarMetBMIAdj$UntarMet_pvalue), method = "BH")
relat_ExprsUntarMetBMIAdj$BMI_FDR <- p.adjust(as.numeric(relat_ExprsUntarMetBMIAdj$BMI_pvalue), method = "BH")
relat_ExprsUntarMetBMIAdj$"UntarMet*BMI_FDR" <- p.adjust(as.numeric(relat_ExprsUntarMetBMIAdj$"UntarMet*BMI_pvalue"), method = "BH")

relat_ProtTarMetBMIAdj <- relat_ProtTarMetBMI
relat_ProtTarMetBMIAdj$TarMet_Q <- qvalue(as.numeric(unlist(relat_ProtTarMetBMIAdj[,"TarMet_pvalue"])))$qvalues
relat_ProtTarMetBMIAdj$BMI_Q <- qvalue(as.numeric(unlist(relat_ProtTarMetBMIAdj[,"BMI_pvalue"])))$qvalues
relat_ProtTarMetBMIAdj$"TarMet*BMI_Q" <- qvalue(as.numeric(unlist(relat_ProtTarMetBMIAdj[,"TarMet*BMI_pvalue"])))$qvalues
relat_ProtTarMetBMIAdj$TarMet_FDR <- p.adjust(as.numeric(relat_ProtTarMetBMIAdj$TarMet_pvalue), method = "BH")
relat_ProtTarMetBMIAdj$BMI_FDR <- p.adjust(as.numeric(relat_ProtTarMetBMIAdj$BMI_pvalue), method = "BH")
relat_ProtTarMetBMIAdj$"TarMet*BMI_FDR" <- p.adjust(as.numeric(relat_ProtTarMetBMIAdj$"TarMet*BMI_pvalue"), method = "BH")

relat_ProtUntarMetBMIAdj <- relat_ProtUntarMetBMI
relat_ProtUntarMetBMIAdj$UntarMet_Q <- qvalue(as.numeric(unlist(relat_ProtUntarMetBMIAdj[,"UntarMet_pvalue"])))$qvalues
relat_ProtUntarMetBMIAdj$BMI_Q <- qvalue(as.numeric(unlist(relat_ProtUntarMetBMIAdj[,"BMI_pvalue"])))$qvalues
relat_ProtUntarMetBMIAdj$"UntarMet*BMI_Q" <- qvalue(as.numeric(unlist(relat_ProtUntarMetBMIAdj[,"UntarMet*BMI_pvalue"])))$qvalues
relat_ProtUntarMetBMIAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_ProtUntarMetBMIAdj$UntarMet_pvalue), method = "BH")
relat_ProtUntarMetBMIAdj$BMI_FDR <- p.adjust(as.numeric(relat_ProtUntarMetBMIAdj$BMI_pvalue), method = "BH")
relat_ProtUntarMetBMIAdj$"UntarMet*BMI_FDR" <- p.adjust(as.numeric(relat_ProtUntarMetBMIAdj$"UntarMet*BMI_pvalue"), method = "BH")
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
relat_ExprsProtHbA1cAdj <- relat_ExprsProtHbA1c
relat_ExprsProtHbA1cAdj$Prot_Q <- qvalue(as.numeric(unlist(relat_ExprsProtHbA1cAdj[,"Prot_pvalue"])))$qvalues
relat_ExprsProtHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(unlist(relat_ExprsProtHbA1cAdj[,"HbA1c_pvalue"])))$qvalues
relat_ExprsProtHbA1cAdj$"Prot*HbA1c_Q" <- qvalue(as.numeric(unlist(relat_ExprsProtHbA1cAdj[,"Prot*HbA1c_pvalue"])))$qvalues
relat_ExprsProtHbA1cAdj$Prot_FDR <- p.adjust(as.numeric(relat_ExprsProtHbA1cAdj$Prot_pvalue), method = "BH")
relat_ExprsProtHbA1cAdj$HbA1c_FDR <- p.adjust(as.numeric(relat_ExprsProtHbA1cAdj$HbA1c_pvalue), method = "BH")
relat_ExprsProtHbA1cAdj$"Prot*HbA1c_FDR" <- p.adjust(as.numeric(relat_ExprsProtHbA1cAdj$"Prot*HbA1c_pvalue"), method = "BH")

relat_ExprsTarMetHbA1cAdj <- relat_ExprsTarMetHbA1c
relat_ExprsTarMetHbA1cAdj$TarMet_Q <- qvalue(as.numeric(unlist(relat_ExprsTarMetHbA1cAdj[,"TarMet_pvalue"])))$qvalues
relat_ExprsTarMetHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(unlist(relat_ExprsTarMetHbA1cAdj[,"HbA1c_pvalue"])))$qvalues
relat_ExprsTarMetHbA1cAdj$"TarMet*HbA1c_Q" <- qvalue(as.numeric(unlist(relat_ExprsTarMetHbA1cAdj[,"TarMet*HbA1c_pvalue"])))$qvalues
relat_ExprsTarMetHbA1cAdj$TarMet_FDR <- p.adjust(as.numeric(relat_ExprsTarMetHbA1cAdj$TarMet_pvalue), method = "BH")
relat_ExprsTarMetHbA1cAdj$HbA1c_FDR <- p.adjust(as.numeric(relat_ExprsTarMetHbA1cAdj$HbA1c_pvalue), method = "BH")
relat_ExprsTarMetHbA1cAdj$"TarMet*HbA1c_FDR" <- p.adjust(as.numeric(relat_ExprsTarMetHbA1cAdj$"TarMet*HbA1c_pvalue"), method = "BH")

relat_ExprsUntarMetHbA1cAdj <- relat_ExprsUntarMetHbA1c
relat_ExprsUntarMetHbA1cAdj$UntarMet_Q <- qvalue(as.numeric(unlist(relat_ExprsUntarMetHbA1cAdj[,"UntarMet_pvalue"])))$qvalues
relat_ExprsUntarMetHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(unlist(relat_ExprsUntarMetHbA1cAdj[,"HbA1c_pvalue"])))$qvalues
relat_ExprsUntarMetHbA1cAdj$"UntarMet*HbA1c_Q" <- qvalue(as.numeric(unlist(relat_ExprsUntarMetHbA1cAdj[,"UntarMet*HbA1c_pvalue"])))$qvalues
relat_ExprsUntarMetHbA1cAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_ExprsUntarMetHbA1cAdj$UntarMet_pvalue), method = "BH")
relat_ExprsUntarMetHbA1cAdj$HbA1c_FDR <- p.adjust(as.numeric(relat_ExprsUntarMetHbA1cAdj$HbA1c_pvalue), method = "BH")
relat_ExprsUntarMetHbA1cAdj$"UntarMet*HbA1c_FDR" <- p.adjust(as.numeric(relat_ExprsUntarMetHbA1cAdj$"UntarMet*HbA1c_pvalue"), method = "BH")

relat_ProtTarMetHbA1cAdj <- relat_ProtTarMetHbA1c
relat_ProtTarMetHbA1cAdj$TarMet_Q <- qvalue(as.numeric(unlist(relat_ProtTarMetHbA1cAdj[,"TarMet_pvalue"])))$qvalues
relat_ProtTarMetHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(unlist(relat_ProtTarMetHbA1cAdj[,"HbA1c_pvalue"])))$qvalues
relat_ProtTarMetHbA1cAdj$"TarMet*HbA1c_Q" <- qvalue(as.numeric(unlist(relat_ProtTarMetHbA1cAdj[,"TarMet*HbA1c_pvalue"])))$qvalues
relat_ProtTarMetHbA1cAdj$TarMet_FDR <- p.adjust(as.numeric(relat_ProtTarMetHbA1cAdj$TarMet_pvalue), method = "BH")
relat_ProtTarMetHbA1cAdj$HbA1c_FDR <- p.adjust(as.numeric(relat_ProtTarMetHbA1cAdj$HbA1c_pvalue), method = "BH")
relat_ProtTarMetHbA1cAdj$"TarMet*HbA1c_FDR" <- p.adjust(as.numeric(relat_ProtTarMetHbA1cAdj$"TarMet*HbA1c_pvalue"), method = "BH")

relat_ProtUntarMetHbA1cAdj <- relat_ProtUntarMetHbA1c
relat_ProtUntarMetHbA1cAdj$UntarMet_Q <- qvalue(as.numeric(unlist(relat_ProtUntarMetHbA1cAdj[,"UntarMet_pvalue"])))$qvalues
relat_ProtUntarMetHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(unlist(relat_ProtUntarMetHbA1cAdj[,"HbA1c_pvalue"])))$qvalues
relat_ProtUntarMetHbA1cAdj$"UntarMet*HbA1c_Q" <- qvalue(as.numeric(unlist(relat_ProtUntarMetHbA1cAdj[,"UntarMet*HbA1c_pvalue"])))$qvalues
relat_ProtUntarMetHbA1cAdj$UntarMet_FDR <- p.adjust(as.numeric(relat_ProtUntarMetHbA1cAdj$UntarMet_pvalue), method = "BH")
relat_ProtUntarMetHbA1cAdj$HbA1c_FDR <- p.adjust(as.numeric(relat_ProtUntarMetHbA1cAdj$HbA1c_pvalue), method = "BH")
relat_ProtUntarMetHbA1cAdj$"UntarMet*HbA1c_FDR" <- p.adjust(as.numeric(relat_ProtUntarMetHbA1cAdj$"UntarMet*HbA1c_pvalue"), method = "BH")

sum(relat_ExprsProt$Prot_FDR < 0.05)
#	1,098,178/6,045,957       18.16%
sum(relat_ExprsTarMet$TarMet_FDR < 0.05)
#	227,634/1,880,244         12.11%
sum(relat_ExprsUntarMet$UntarMet_FDR < 0.05)
#	227,909/3,776,697         6.03%
# Exprs ~ Metab (merged): 455,543/5,656,941	8.05%
sum(relat_ProtTarMet$TarMet_FDR < 0.05)
#	20,056/43,268             46.35%
sum(relat_ProtUntarMet$UntarMet_FDR < 0.05)
#	2/86,909                  0%
#Prot ~ Metab (merged): 20,058/130,177	15.41%

#	TOTAL PAIRS: 11,833,075
#	SIGNIFICANT PAIRS: 1,941,798  16.41%
#	SIGNIFICANT PAIRS (FDR): 1,573,779    13.3%
# TOTAL MOLECULES: 16,931
# UNIQUE MOLECULES IN SIGNIF PAIRS: 16,931 100%
# UNIQUE MOLECULES IN SIGNIF PAIRS (FDR): 16,927 99.9%
#------------------
#     1. SEX
#------------------
sum(relat_ExprsProtSexAdj$"Prot*Sex_FDR" < 0.05)
#  3,243,676 ( associations)
sum(relat_ExprsTarMetSexAdj$"TarMet*Sex_FDR" < 0.05)
#	1,021,706 ( associations)
sum(relat_ExprsUntarMetSexAdj$"UntarMet*Sex_FDR"< 0.05)
#	2,070,835 ( associations)
sum(relat_ProtTarMetSexAdj$"TarMet*Sex_FDR" < 0.05)
#	28,642 ( associations)
sum(relat_ProtUntarMetSexAdj$"UntarMet*Sex_FDR" < 0.05)
#	4,575 ( associations)

#------------------
#     2. BMI
#------------------
sum(relat_ExprsProtBMIAdj$"Prot*BMI_FDR" < 0.05)
# 3,246,257 ( associations)
sum(relat_ExprsTarMetBMIAdj$"TarMet*BMI_FDR" < 0.05)
# 1,022,594 ( associations)
sum(relat_ExprsUntarMetBMIAdj$"UntarMet*BMI_FDR"< 0.05)
# 2,070,374 ( associations)
sum(relat_ProtTarMetBMIAdj$"TarMet*BMI_FDR" < 0.05)
# 28,540 ( associations)
sum(relat_ProtUntarMetBMIAdj$"UntarMet*BMI_FDR" < 0.05)
# 4,535 ( associations)

#------------------
#     3. HbA1c
#------------------
sum(relat_ExprsProtHbA1cAdj$"Prot*HbA1c_FDR" < 0.05)
#	3,245,875 ( associations)
sum(relat_ExprsTarMetHbA1cAdj$"TarMet*HbA1c_FDR" < 0.05)
#	1,023,938 ( associations)
sum(relat_ExprsUntarMetHbA1cAdj$"UntarMet*HbA1c_FDR"< 0.05)
#	2,070,475 ( associations)
sum(relat_ProtTarMetHbA1cAdj$"TarMet*HbA1c_FDR" < 0.05)
#	28,526 ( associations)
sum(relat_ProtUntarMetHbA1cAdj$"UntarMet*HbA1c_FDR" < 0.05)
#	4,593 ( associations)

#------------------
#     4. AGE
#------------------
sum(relat_ExprsProtAgeAdj$"Prot*Age_FDR" < 0.05)
# 2,765,270 ( associations)
sum(relat_ExprsTarMetAgeAdj$"TarMet*Age_FDR" < 0.05)
# 857,940 ( associations)
sum(relat_ExprsUntarMetAgeAdj$"UntarMet*Age_FDR"< 0.05)
# 1,718,187 ( associations)
sum(relat_ProtTarMetAgeAdj$"TarMet*Age_FDR" < 0.05)
# 24,173 ( associations)
sum(relat_ProtUntarMetAgeAdj$"UntarMet*Age_FDR" < 0.05)
# 640 ( associations)
#------------------------------------------------------------------------------------------------------------------------------------------------


##		BARPLOT
pairInter_counts <- matrix(0,ncol=4,nrow=6)
colnames(pairInter_counts) <- c("Age","Sex","BMI","HbA1c")
rownames(pairInter_counts) <- c("Gene-Protein Interaction","Gene-Protein Interaction & Association","Gene-Metabolite Interaction","Gene-Metabolite Interaction & Association","Protein-Metabolite Interaction","Protein-Metabolite Interaction & Association")
pairInter_counts[1,1] <- sum(relat_ExprsProtAgeAdj$"Prot*Age_FDR"<0.05)
pairInter_counts[2,1] <- sum(relat_ExprsProtAgeAdj$"Prot*Age_FDR"<0.05&relat_ExprsProtAgeAdj$"Prot_FDR"<0.05)
pairInter_counts[3,1] <- sum(relat_ExprsTarMetAgeAdj$"TarMet*Age_FDR"<0.05)+sum(relat_ExprsUntarMetAgeAdj$"UntarMet*Age_FDR"<0.05)
pairInter_counts[4,1] <- sum(relat_ExprsTarMetAgeAdj$"TarMet*Age_FDR"<0.05&relat_ExprsTarMetAgeAdj$"TarMet_FDR"<0.05)+sum(relat_ExprsUntarMetAgeAdj$"UntarMet*Age_FDR"<0.05&relat_ExprsUntarMetAgeAdj$"UntarMet_FDR"<0.05)
pairInter_counts[5,1] <- sum(relat_ProtTarMetAgeAdj$"TarMet*Age_FDR"<0.05)+sum(relat_ProtUntarMetAgeAdj$"UntarMet*Age_FDR"<0.05)
pairInter_counts[6,1] <- sum(relat_ProtTarMetAgeAdj$"TarMet*Age_FDR"<0.05&relat_ProtTarMetAgeAdj$"TarMet_FDR"<0.05)+sum(relat_ProtUntarMetAgeAdj$"UntarMet*Age_FDR"<0.05&relat_ProtUntarMetAgeAdj$"UntarMet_FDR"<0.05)
pairInter_counts[1,2] <- sum(relat_ExprsProtSexAdj$"Prot*Sex_FDR"<0.05)
pairInter_counts[2,2] <- sum(relat_ExprsProtSexAdj$"Prot*Sex_FDR"<0.05&relat_ExprsProtSexAdj$"Prot_FDR"<0.05)
pairInter_counts[3,2] <- sum(relat_ExprsTarMetSexAdj$"TarMet*Sex_FDR"<0.05)+sum(relat_ExprsUntarMetSexAdj$"UntarMet*Sex_FDR"<0.05)
pairInter_counts[4,2] <- sum(relat_ExprsTarMetSexAdj$"TarMet*Sex_FDR"<0.05&relat_ExprsTarMetSexAdj$"TarMet_FDR"<0.05)+sum(relat_ExprsUntarMetSexAdj$"UntarMet*Sex_FDR"<0.05&relat_ExprsUntarMetSexAdj$"UntarMet_FDR"<0.05)
pairInter_counts[5,2] <- sum(relat_ProtTarMetSexAdj$"TarMet*Sex_FDR"<0.05)+sum(relat_ProtUntarMetSexAdj$"UntarMet*Sex_FDR"<0.05)
pairInter_counts[6,2] <- sum(relat_ProtTarMetSexAdj$"TarMet*Sex_FDR"<0.05&relat_ProtTarMetSexAdj$"TarMet_FDR"<0.05)+sum(relat_ProtUntarMetSexAdj$"UntarMet*Sex_FDR"<0.05&relat_ProtUntarMetSexAdj$"UntarMet_FDR"<0.05)
pairInter_counts[1,3] <- sum(relat_ExprsProtBMIAdj$"Prot*BMI_FDR"<0.05)
pairInter_counts[2,3] <- sum(relat_ExprsProtBMIAdj$"Prot*BMI_FDR"<0.05&relat_ExprsProtBMIAdj$"Prot_FDR"<0.05)
pairInter_counts[3,3] <- sum(relat_ExprsTarMetBMIAdj$"TarMet*BMI_FDR"<0.05)+sum(relat_ExprsUntarMetBMIAdj$"UntarMet*BMI_FDR"<0.05)
pairInter_counts[4,3] <- sum(relat_ExprsTarMetBMIAdj$"TarMet*BMI_FDR"<0.05&relat_ExprsTarMetBMIAdj$"TarMet_FDR"<0.05)+sum(relat_ExprsUntarMetBMIAdj$"UntarMet*BMI_FDR"<0.05&relat_ExprsUntarMetBMIAdj$"UntarMet_FDR"<0.05)
pairInter_counts[5,3] <- sum(relat_ProtTarMetBMIAdj$"TarMet*BMI_FDR"<0.05)+sum(relat_ProtUntarMetBMIAdj$"UntarMet*BMI_FDR"<0.05)
pairInter_counts[6,3] <- sum(relat_ProtTarMetBMIAdj$"TarMet*BMI_FDR"<0.05&relat_ProtTarMetBMIAdj$"TarMet_FDR"<0.05)+sum(relat_ProtUntarMetBMIAdj$"UntarMet*BMI_FDR"<0.05&relat_ProtUntarMetBMIAdj$"UntarMet_FDR"<0.05)
pairInter_counts[1,4] <- sum(relat_ExprsProtHbA1cAdj$"Prot*HbA1c_FDR"<0.05)
pairInter_counts[2,4] <- sum(relat_ExprsProtHbA1cAdj$"Prot*HbA1c_FDR"<0.05&relat_ExprsProtHbA1cAdj$"Prot_FDR"<0.05)
pairInter_counts[3,4] <- sum(relat_ExprsTarMetHbA1cAdj$"TarMet*HbA1c_FDR"<0.05)+sum(relat_ExprsUntarMetHbA1cAdj$"UntarMet*HbA1c_FDR"<0.05)
pairInter_counts[4,4] <- sum(relat_ExprsTarMetHbA1cAdj$"TarMet*HbA1c_FDR"<0.05&relat_ExprsTarMetHbA1cAdj$"TarMet_FDR"<0.05)+sum(relat_ExprsUntarMetHbA1cAdj$"UntarMet*HbA1c_FDR"<0.05&relat_ExprsUntarMetHbA1cAdj$"UntarMet_FDR"<0.05)
pairInter_counts[5,4] <- sum(relat_ProtTarMetHbA1cAdj$"TarMet*HbA1c_FDR"<0.05)+sum(relat_ProtUntarMetHbA1cAdj$"UntarMet*HbA1c_FDR"<0.05)
pairInter_counts[6,4] <- sum(relat_ProtTarMetHbA1cAdj$"TarMet*HbA1c_FDR"<0.05&relat_ProtTarMetHbA1cAdj$"TarMet_FDR"<0.05)+sum(relat_ProtUntarMetHbA1cAdj$"UntarMet*HbA1c_FDR"<0.05&relat_ProtUntarMetHbA1cAdj$"UntarMet_FDR"<0.05)

pairInter_percentages <- matrix(0,ncol=4,nrow=6)
colnames(pairInter_percentages) <- c("Age","Sex","BMI","HbA1c")
rownames(pairInter_percentages) <- c("Gene-Protein Interaction","Gene-Protein Interaction & Association","Gene-Metabolite Interaction","Gene-Metabolite Interaction & Association","Protein-Metabolite Interaction","Protein-Metabolite Interaction & Association")
pairInter_percentages[1,1] <- pairInter_counts[1,1]*100/6045957
pairInter_percentages[1,2] <- pairInter_counts[1,2]*100/6045957
pairInter_percentages[1,3] <- pairInter_counts[1,3]*100/6045957
pairInter_percentages[1,4] <- pairInter_counts[1,4]*100/6045957
pairInter_percentages[2,1] <- pairInter_counts[2,1]*100/6045957
pairInter_percentages[2,2] <- pairInter_counts[2,2]*100/6045957
pairInter_percentages[2,3] <- pairInter_counts[2,3]*100/6045957
pairInter_percentages[2,4] <- pairInter_counts[2,4]*100/6045957
pairInter_percentages[3,1] <- pairInter_counts[3,1]*100/5656941
pairInter_percentages[3,2] <- pairInter_counts[3,2]*100/5656941
pairInter_percentages[3,3] <- pairInter_counts[3,3]*100/5656941
pairInter_percentages[3,4] <- pairInter_counts[3,4]*100/5656941
pairInter_percentages[4,1] <- pairInter_counts[4,1]*100/5656941
pairInter_percentages[4,2] <- pairInter_counts[4,2]*100/5656941
pairInter_percentages[4,3] <- pairInter_counts[4,3]*100/5656941
pairInter_percentages[4,4] <- pairInter_counts[4,4]*100/5656941
pairInter_percentages[5,1] <- pairInter_counts[5,1]*100/130177
pairInter_percentages[5,2] <- pairInter_counts[5,2]*100/130177
pairInter_percentages[5,3] <- pairInter_counts[5,3]*100/130177
pairInter_percentages[5,4] <- pairInter_counts[5,4]*100/130177
pairInter_percentages[6,1] <- pairInter_counts[6,1]*100/130177
pairInter_percentages[6,2] <- pairInter_counts[6,2]*100/130177
pairInter_percentages[6,3] <- pairInter_counts[6,3]*100/130177
pairInter_percentages[6,4] <- pairInter_counts[6,4]*100/130177

pairInter_percentages2 <- data.frame(Interactions = c(pairInter_percentages[1:6,1],pairInter_percentages[1:6,2],pairInter_percentages[1:6,3],pairInter_percentages[1:6,4]), 
                                     Covariate = rep(c("Age", "Sex", "BMI","HbA1c"), each = 6),
                                     Pair = c("Gene-Protein Interaction","Gene-Protein Interaction & Association","Gene-Metabolite Interaction","Gene-Metabolite Interaction & Association","Protein-Metabolite Interaction","Protein-Metabolite Interaction & Association"))
pairInter_percentages2$Group <- rep(c(rep("Gene-Protein",2),rep("Gene-Metabolite",2),rep("Protein-Metabolite",2)),4)
#pairInter_percentages2$Set <- rep(c("all","sub"),12)
ggplot(pairInter_percentages2, aes(x = Group, y =Interactions, fill = factor(Pair, levels = c("Gene-Protein Interaction","Gene-Protein Interaction & Association","Gene-Metabolite Interaction","Gene-Metabolite Interaction & Association","Protein-Metabolite Interaction","Protein-Metabolite Interaction & Association")))) +
  geom_bar(stat = "identity", position = "identity", width = 0.7) +
  #geom_bar(position="stack", stat="identity") +
  labs(x=NULL, y = "% Significant Interactions", fill = NULL) +
  scale_fill_manual(values=c("darkmagenta","mediumorchid3","darkgoldenrod","goldenrod2","lightsteelblue4","lightsteelblue3")) +
  facet_wrap(~Covariate) +
  theme(strip.text.x = element_text(size = 12, face = "bold" ),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14, angle = 90),
        axis.title = element_text(size = 18),
        legend.text = element_text(size = 16))


##	UpSETR PLOTS
G_P_Plot <- data.frame(paste(relat_ExprsProtAgeAdj$GeneID,relat_ExprsProtAgeAdj$ProteinID,sep="_"),relat_ExprsProtAgeAdj$"Prot*Age_FDR",relat_ExprsProtSexAdj$"Prot*Sex_FDR",relat_ExprsProtBMIAdj$"Prot*BMI_FDR",relat_ExprsProtHbA1cAdj$"Prot*HbA1c_FDR")
colnames(G_P_Plot) <- c("PAIR","Prot*Age_FDR","Prot*BMI_FDR","Prot*Sex_FDR","Prot*HbA1c_FDR")
G_P_Plot$Age <- 0
G_P_Plot$Age[G_P_Plot$"Prot*Age_FDR" < 0.05] <- 1
G_P_Plot$BMI <- 0
G_P_Plot$BMI[G_P_Plot$"Prot*BMI_FDR" < 0.05] <- 1
G_P_Plot$Sex <- 0
G_P_Plot$Sex[G_P_Plot$"Prot*Sex_FDR" < 0.05] <- 1
G_P_Plot$HbA1c <- 0
G_P_Plot$HbA1c[G_P_Plot$"Prot*HbA1c_FDR" < 0.05] <- 1
G_P_Plot[,c(2:5)] <- NULL
upset(G_P_Plot, sets = c("Age", "BMI", "Sex", "HbA1c"), main.bar.color="darkmagenta" , matrix.color="black", sets.bar.color = "grey",
      keep.order = TRUE, mainbar.y.label = "Common Molecule Pairs", sets.x.label = "Significant Interaction Count", 
      text.scale = c(1.7,1.3,1.3,1.3,1.4,1.5), order.by = "degree", empty.intersections = "on")
#------------------------------------------------------------------------------------------------------
gene_metab <- c(paste(relat_ExprsTarMetAgeAdj$GeneID,relat_ExprsTarMetAgeAdj$MetaboliteID,sep="_"),paste(relat_ExprsUntarMetAgeAdj$GeneID,relat_ExprsUntarMetAgeAdj$BIOCHEMICAL2,sep="_"))
gene_metab_Age <- c(relat_ExprsTarMetAgeAdj$"TarMet*Age_FDR",relat_ExprsUntarMetAgeAdj$"UntarMet*Age_FDR")
gene_metab_Sex <- c(relat_ExprsTarMetSexAdj$"TarMet*Sex_FDR",relat_ExprsUntarMetSexAdj$"UntarMet*Sex_FDR")
gene_metab_BMI <- c(relat_ExprsTarMetBMIAdj$"TarMet*BMI_FDR",relat_ExprsUntarMetBMIAdj$"UntarMet*BMI_FDR")
gene_metab_HbA1c <- c(relat_ExprsTarMetHbA1cAdj$"TarMet*HbA1c_FDR",relat_ExprsUntarMetHbA1cAdj$"UntarMet*HbA1c_FDR")

G_M_Plot <- data.frame(gene_metab,gene_metab_Age,gene_metab_Sex,gene_metab_BMI,gene_metab_HbA1c)
colnames(G_M_Plot) <- c("PAIR","Met*Age_FDR","Met*BMI_FDR","Met*Sex_FDR","Met*HbA1c_FDR")
G_M_Plot$Age <- 0
G_M_Plot$Age[G_M_Plot$"Met*Age_FDR" < 0.05] <- 1
G_M_Plot$BMI <- 0
G_M_Plot$BMI[G_M_Plot$"Met*BMI_FDR" < 0.05] <- 1
G_M_Plot$Sex <- 0
G_M_Plot$Sex[G_M_Plot$"Met*Sex_FDR" < 0.05] <- 1
G_M_Plot$HbA1c <- 0
G_M_Plot$HbA1c[G_M_Plot$"Met*HbA1c_FDR" < 0.05] <- 1
G_M_Plot[,c(2:5)] <- NULL
upset(G_M_Plot, sets = c("Age", "BMI", "Sex", "HbA1c"), main.bar.color="darkgoldenrod" , matrix.color="black", sets.bar.color = "grey",
      keep.order = TRUE, mainbar.y.label = "Common Molecule Pairs", sets.x.label = "Significant Interaction Count", 
      text.scale = c(1.7,1.3,1.3,1.3,1.4,1.5), order.by = "degree", empty.intersections = "on")
#------------------------------------------------------------------------------------------------------
prot_metab <- c(paste(relat_ProtTarMetAgeAdj$ProteinID,relat_ProtTarMetAgeAdj$MetaboliteID,sep="_"),paste(relat_ProtUntarMetAgeAdj$GeneID,relat_ProtUntarMetAgeAdj$BIOCHEMICAL2,sep="_"))
prot_metab_Age <- c(relat_ProtTarMetAgeAdj$"TarMet*Age_FDR",relat_ProtUntarMetAgeAdj$"UntarMet*Age_FDR")
prot_metab_Sex <- c(relat_ProtTarMetSexAdj$"TarMet*Sex_FDR",relat_ProtUntarMetSexAdj$"UntarMet*Sex_FDR")
prot_metab_BMI <- c(relat_ProtTarMetBMIAdj$"TarMet*BMI_FDR",relat_ProtUntarMetBMIAdj$"UntarMet*BMI_FDR")
prot_metab_HbA1c <- c(relat_ProtTarMetHbA1cAdj$"TarMet*HbA1c_FDR",relat_ProtUntarMetHbA1cAdj$"UntarMet*HbA1c_FDR")

P_M_Plot <- data.frame(prot_metab,prot_metab_Age,prot_metab_Sex,prot_metab_BMI,prot_metab_HbA1c)
colnames(P_M_Plot) <- c("PAIR","Met*Age_FDR","Met*BMI_FDR","Met*Sex_FDR","Met*HbA1c_FDR")
P_M_Plot$Age <- 0
P_M_Plot$Age[P_M_Plot$"Met*Age_FDR" < 0.05] <- 1
P_M_Plot$BMI <- 0
P_M_Plot$BMI[P_M_Plot$"Met*BMI_FDR" < 0.05] <- 1
P_M_Plot$Sex <- 0
P_M_Plot$Sex[P_M_Plot$"Met*Sex_FDR" < 0.05] <- 1
P_M_Plot$HbA1c <- 0
P_M_Plot$HbA1c[P_M_Plot$"Met*HbA1c_FDR" < 0.05] <- 1
P_M_Plot[,c(2:5)] <- NULL
upset(P_M_Plot, sets = c("Age", "BMI", "Sex", "HbA1c"), main.bar.color="lightsteelblue4" , matrix.color="black", sets.bar.color = "grey",
      keep.order = TRUE, mainbar.y.label = "Common Molecule Pairs", sets.x.label = "Significant Interaction Count", 
      text.scale = c(1.7,1.3,1.3,1.3,1.4,1.5), order.by = "degree", empty.intersections = "on")
	  
	  
##	INTERACTION PLOTS
panel <- rbind(c('Olink_CARDIOMETABOLIC','C-Met','CAM'), c('Olink_CARDIOVASCULAR_II','CVD_II','CVD2'),
               c('Olink_CARDIOVASCULAR_III','CVD_III','CVD3'), c('Olink_DEVELOPMENT','Dev','DEV'),
               c('Olink_METABOLISM','Met','MET'))
raw_genesE2 <- raw_genesE[,which(colnames(raw_genesE) %in% colnames(raw_proteins_Imp))]
covar_Exprsf2 <- covar_Exprsf[which(colnames(raw_genesE) %in% colnames(raw_proteins_Imp)),]
raw_met_targE2 <- raw_met_targE[,which(colnames(raw_met_targE) %in% colnames(raw_proteins_Imp))]
covar_Targf2 <- covar_Targf[which(colnames(raw_met_targE) %in% colnames(raw_proteins_Imp)),]

for(i in 1:nrow(covar_Exprsf2)){
  if(covar_Exprsf2$BMI[i] <= 18.5){covar_Exprsf2$BMI_cat[i] <- "Underweight"}
  else if(covar_Exprsf2$BMI[i] > 18.5 & covar_Exprsf2$BMI[i] < 25){covar_Exprsf2$BMI_cat[i] <- "Normal"}
  else if(covar_Exprsf2$BMI[i] >= 25 & covar_Exprsf2$BMI[i] < 30){covar_Exprsf2$BMI_cat[i] <- "Overweight"}
  else{covar_Exprsf2$BMI_cat[i] <- "Obese"}}

for(i in 1:nrow(covar_Exprsf2)){
  if(covar_Exprsf2$HbA1c[i] < 39.7){covar_Exprsf2$HbA1c_cat[i] <- "Low"}
  else{covar_Exprsf2$HbA1c_cat[i] <- "High"}}

table(covar_Exprsf2$BMI_cat)
# Normal       Obese  Overweight Underweight 
#   568          963       1489          3
table(covar_Exprsf2$HbA1c_cat)
# High  Low
# 1184  1839

#	20% of 3,023 is 605 samples
rownames(covar_Exprsf2) <- colnames(raw_genesE2)
SampleIDs_covar_Exprsf2_HbA1cLOW <- cbind(rownames(covar_Exprsf2[order(covar_Exprsf2[,"HbA1c"]),])[1:605],covar_Exprsf2[order(covar_Exprsf2[,"HbA1c"]),"HbA1c"][1:605])
SampleIDs_covar_Exprsf2_HbA1cHIGH <- cbind(rownames(covar_Exprsf2[order(covar_Exprsf2[,"HbA1c"],decreasing=TRUE),])[1:605],covar_Exprsf2[order(covar_Exprsf2[,"HbA1c"],decreasing=TRUE),"HbA1c"][1:605])
covar_Exprsf2$ID <- rownames(covar_Exprsf2)
covar_Exprsf2_HbA1cSub <- subset(covar_Exprsf2, subset=((ID %in% SampleIDs_covar_Exprsf2_HbA1cLOW[,1])|(ID %in% SampleIDs_covar_Exprsf2_HbA1cHIGH[,1])))
covar_Exprsf2_HbA1cSub$HbA1c_group <- ifelse(covar_Exprsf2_HbA1cSub$HbA1c > 39,"High","Low")
SampleIDs_Prot_TOP20pc <- covar_Exprsf2_HbA1cSub$ID
raw_genesE2_HbA1cSub <- raw_genesE2[,SampleIDs_Prot_TOP20pc]
raw_proteins_Imp_HbA1cSub <- raw_proteins_Imp[,SampleIDs_Prot_TOP20pc]
covar_Protf$ID <- colnames(raw_proteins_Imp)
covar_Protf_HbA1cSub <- subset(covar_Protf, ID %in% SampleIDs_Prot_TOP20pc)
covar_Targf2$ID <- colnames(raw_met_targE2)
raw_met_targE2_HbA1cSub <- raw_met_targE2[,SampleIDs_Prot_TOP20pc]
#---------------------------------------------------------------------------
SampleIDs_covar_Exprsf2_AgeLOW <- cbind(rownames(covar_Exprsf2[order(covar_Exprsf2[,"Age"]),])[1:605],covar_Exprsf2[order(covar_Exprsf2[,"Age"]),"Age"][1:605])
SampleIDs_covar_Exprsf2_AgeHIGH <- cbind(rownames(covar_Exprsf2[order(covar_Exprsf2[,"Age"],decreasing=TRUE),])[1:605],covar_Exprsf2[order(covar_Exprsf2[,"Age"],decreasing=TRUE),"Age"][1:605])
covar_Exprsf2_AgeSub <- subset(covar_Exprsf2, subset=((ID %in% SampleIDs_covar_Exprsf2_AgeLOW[,1])|(ID %in% SampleIDs_covar_Exprsf2_AgeHIGH[,1])))
covar_Exprsf2_AgeSub$Age_group <- ifelse(covar_Exprsf2_AgeSub$Age > 61,"Older","Younger")
SampleIDs_AGE_TOP20pc <- covar_Exprsf2_AgeSub$ID
raw_genesE2_AgeSub <- raw_genesE2[,SampleIDs_AGE_TOP20pc]
raw_proteins_Imp_AgeSub <- raw_proteins_Imp[,SampleIDs_AGE_TOP20pc]
covar_Protf_AgeSub <- subset(covar_Protf, ID %in% SampleIDs_AGE_TOP20pc)
raw_met_targE2_AgeSub <- raw_met_targE2[,SampleIDs_AGE_TOP20pc]
covar_Targf2_AgeSub <- subset(covar_Targf2, ID %in% SampleIDs_AGE_TOP20pc)
#-------------------------------------------------------------------------------------------------------------
colnames(relat_ExprsProtSexAdj)[4] <- "ProteinName"
relat_ExprsProtSexAdj2 <- relat_ExprsProtSexAdj[which(relat_ExprsProtSexAdj$"Prot*Sex_FDR" < 0.05),]
dim(subset(relat_ExprsProtSexAdj2, subset=(ProteinName=="IGFBP1")))
# [1] 8624   19
IGFBP1_genesS <- unique(subset(relat_ExprsProtSexAdj2, subset=(ProteinName=="IGFBP1"))[,2])
"NOC2L" %in% IGFBP1_genesS$GeneName
"RAB1A" %in% unique(subset(relat_ExprsProtSexAdj2, subset=(ProteinName=="GDF15"))[,2])$GeneName

colnames(relat_ProtTarMetSexAdj)[2] <- "ProteinName"
relat_ProtTarMetSexAdj2 <- relat_ProtTarMetSexAdj[which(relat_ProtTarMetSexAdj$"TarMet*Sex_FDR" < 0.05),]
dim(subset(relat_ProtTarMetSexAdj2, subset=(ProteinName=="GDF15")))
# [1] 116   19
IGFBP2_metabS <- unique(subset(relat_ProtTarMetSexAdj2, subset=(ProteinName=="IGFBP2"))[,3])
"C0" %in% IGFBP2_metabS$MetaboliteID
"C0" %in% unique(subset(relat_ProtTarMetSexAdj2, subset=(ProteinName=="IGFBP1"))[,3])

relat_ExprsTarMetSexAdj2 <- relat_ExprsTarMetSexAdj[which(relat_ExprsTarMetSexAdj$"TarMet*Sex_FDR" < 0.05),]
"C0" %in% unique(subset(relat_ExprsTarMetSexAdj2, subset=(GeneName=="FAM63A"))[,3])

colnames(relat_ExprsProtBMIAdj)[4] <- "ProteinName"
relat_ExprsProtBMIAdj2 <- relat_ExprsProtBMIAdj[which(relat_ExprsProtBMIAdj$"Prot*BMI_FDR" < 0.05),]
dim(subset(relat_ExprsProtBMIAdj2, subset=(ProteinName=="IGFBP1")))
# [1] 8594   19
IGFBP1_genesB <- unique(subset(relat_ExprsProtBMIAdj2, subset=(ProteinName=="IGFBP1"))[,2])
"POGZ" %in% IGFBP1_genesB$GeneName

colnames(relat_ExprsProtHbA1cAdj)[4] <- "ProteinName"
relat_ExprsProtHbA1cAdj2 <- relat_ExprsProtHbA1cAdj[which(relat_ExprsProtHbA1cAdj$"Prot*HbA1c_FDR" < 0.05),]
dim(subset(relat_ExprsProtHbA1cAdj2, subset=(ProteinName=="IL6")))
# [1] 7719   19
IL6_genesB <- unique(subset(relat_ExprsProtHbA1cAdj2, subset=(ProteinName=="IL6"))[,2])
"FAM63A" %in% IL6_genesB$GeneName
"IGFBP2" %in% IL6_genesB$GeneName

relat_ExprsTarMetBMIAdj2 <- relat_ExprsTarMetBMIAdj[which(relat_ExprsTarMetBMIAdj$"TarMet*BMI_FDR" < 0.05),]
"C0" %in% unique(subset(relat_ExprsTarMetBMIAdj2, subset=(GeneName=="FAM63A"))[,3])

relat_ExprsTarMetHbA1cAdj2 <- relat_ExprsTarMetHbA1cAdj[which(relat_ExprsTarMetHbA1cAdj$"TarMet*HbA1c_FDR" < 0.05),]
"C0" %in% unique(subset(relat_ExprsTarMetHbA1cAdj2, subset=(GeneName=="FAM63A"))[,3])

colnames(relat_ProtTarMetBMIAdj)[2] <- "ProteinName"
relat_ProtTarMetBMIAdj2 <- relat_ProtTarMetBMIAdj[which(relat_ProtTarMetBMIAdj$"TarMet*BMI_FDR" < 0.05),]
"C0" %in% unique(subset(relat_ProtTarMetBMIAdj2, subset=(ProteinName=="GDF15"))[,3])$MetaboliteID

relat_ProtTarMetHbA1cAdj2 <- relat_ProtTarMetHbA1cAdj[which(relat_ProtTarMetHbA1cAdj$"TarMet*HbA1c_FDR" < 0.05),]
"C0" %in% unique(subset(relat_ProtTarMetHbA1cAdj2, subset=(GeneName=="IGFBP2"))[,3])$MetaboliteID

relat_ProtTarMetAgeAdj2 <- relat_ProtTarMetAgeAdj[which(relat_ProtTarMetAgeAdj$"TarMet*Age_FDR" < 0.05),]
"C0" %in% unique(subset(relat_ProtTarMetAgeAdj2, subset=(GeneName=="GDF15"))[,3])$MetaboliteID

relat_ExprsTarMetAgeAdj2 <- relat_ExprsTarMetAgeAdj[which(relat_ExprsTarMetAgeAdj$"TarMet*Age_FDR" < 0.05),]
"C0" %in% unique(subset(relat_ExprsTarMetAgeAdj2, subset=(GeneName=="FAM63A"))[,3])

colnames(relat_ExprsProtAgeAdj)[4] <- "ProteinName"
relat_ExprsProtAgeAdj2 <- relat_ExprsProtAgeAdj[which(relat_ExprsProtAgeAdj$"Prot*Age_FDR" < 0.05),]
"METRNL" %in% unique(subset(relat_ExprsProtAgeAdj2, subset=(GeneName =="RAB1A"))[,4])$ProteinName

colnames(relat_ExprsProtBMIAdj)[4] <- "ProteinName"
relat_ExprsProtBMIAdj2 <- relat_ExprsProtBMIAdj[which(relat_ExprsProtBMIAdj$"Prot*BMI_FDR" < 0.05),]
"RAB1A" %in% unique(subset(relat_ExprsProtBMIAdj2, subset=(ProteinName=="IGFBP2"))[,2])$GeneName
"FAM63A" %in% unique(subset(relat_ExprsProtBMIAdj2, subset=(ProteinName=="METRNL"))[,2])$GeneName

colnames(relat_ExprsProtHbA1cAdj)[4] <- "ProteinName"
relat_ExprsProtHbA1cAdj2 <- relat_ExprsProtHbA1cAdj[which(relat_ExprsProtHbA1cAdj$"Prot*HbA1c_FDR" < 0.05),]
dim(subset(relat_ExprsProtHbA1cAdj2, subset=(ProteinName=="IL6")))
# [1] 7719   19
IL6_genesB <- unique(subset(relat_ExprsProtHbA1cAdj2, subset=(ProteinName=="IL6"))[,2])
"FAM63A" %in% IL6_genesB$GeneName
#-------------------------------------------------------------------------------------------------------------
which(anno_genes$GeneName == " RAB1A")
# 1974
which(anno_genes$GeneName == " RFX5")
# 1012
which(anno_genes$GeneName == " PRUNE")
# 995
which(anno_genes$GeneName == " TSPAN3")
# 11257
which(anno_genes$GeneName == " POGZ")
# 1015
which(anno_genes$GeneName == " NOC2L")
# 12
which(anno_genes$GeneName == " IGFBP2")
# 2571
which(anno_genes$GeneName == " FAM63A")
# 994
which(anno_proteins$GeneName == " IGFBP1")
# 177
which(anno_proteins$GeneName == " IGFBP2")
# 221
which(anno_proteins$GeneName == " GDF15")
# 172
which(anno_proteins$GeneName == " IL6")
# 79
which(anno_proteins$GeneName == " APOM")
# 55
which(anno_proteins$GeneName == " METRNL")
# 370
which(anno_met_targ$Metabolite == "C0") # L-Carnitine
# 1

plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[177,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(as.double(raw_genesE2[12,])) ~ 1 + . , data = covar_Exprsf2[c(1,3:5,7:11)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp[177,])) ~ 1 + . + plate, data = covar_Protf[c(1,3:5,7)])), Sex =covar_Exprsf2$Sex)
pdf(file = "NOC2L_IGFBP1_Sex.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe2, y = MolPhe1, color  = Sex)) + geom_point(alpha=0.4,shape=20) + scale_color_manual(values=c("darkred", "darkslateblue")) +
  stat_smooth(data=subset(data, Sex == " female"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkred") +
  stat_smooth(data=subset(data, Sex == " male"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslateblue") +
  labs(title="NOC2L gene expression VS IGFBP1 levels", x ="IGFBP1", y = "NOC2L") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[172,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(as.double(raw_genesE2[1974,])) ~ 1 + . , data = covar_Exprsf2[c(1,3:5,7:11)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp[172,])) ~ 1 + . + plate, data = covar_Protf[c(1,3:5,7)])), Sex =covar_Exprsf2$Sex)
pdf(file = "RAB1A_GDF15_Sex.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe2, y = MolPhe1, color  = Sex)) + geom_point(alpha=0.4,shape=20) + scale_color_manual(values=c("darkred", "darkslateblue")) +
  stat_smooth(data=subset(data, Sex == " female"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkred") +
  stat_smooth(data=subset(data, Sex == " male"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslateblue") +
  labs(title="RAB1A gene expression VS GDF15 levels", x ="GDF15", y = "RAB1A") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[221,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(t(raw_met_targE2[1,])) ~ 1 + . , data = covar_Targf2[c(1,3:5,7:9)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp[221,])) ~ 1 + . + plate, data = covar_Protf[c(1,3:5,7)])), Sex =covar_Exprsf2$Sex)
pdf(file = "IGFBP2_Carnitine_Sex.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe1, y = MolPhe2, color  = Sex)) + geom_point(alpha=0.4,shape=20) + scale_color_manual(values=c("darkred", "darkslateblue")) +
  stat_smooth(data=subset(data, Sex == " female"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkred") +
  stat_smooth(data=subset(data, Sex == " male"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslateblue") +
  labs(title="IGFBP2 levels VS L-Carnitine levels", x ="L-Carnitine", y = "IGFBP2") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[177,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(as.double(raw_genesE2[1015,])) ~ 1 + . , data = covar_Exprsf2[c(1:2,4:5,7:11)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp[177,])) ~ 1 + . + plate, data = covar_Protf[c(1:2,4:5,7)])), BMI =covar_Exprsf2$BMI_cat)
pdf(file = "POGZ_IGFBP1_BMI.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe2, y = MolPhe1, color  = BMI)) + geom_point(alpha=0.5,shape=20) +
  scale_color_manual(name = "BMI", values=setNames(c("darkslateblue", "darkslategray", "goldenrod4", "maroon4"),c("Underweight", "Normal", "Overweight", "Obese"))) +
  stat_smooth(data=subset(data, BMI == "Normal"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslategray") +
  #stat_smooth(data=subset(data, BMI == "Overweight"), method = "lm", se = FALSE, fullrange = TRUE, col = "goldenrod4") +
  stat_smooth(data=subset(data, BMI == "Obese"), method = "lm", se = FALSE, fullrange = TRUE, col = "maroon4") +
  labs(title="POGZ gene expression VS IGFBP1 levels", x ="IGFBP1", y = "POGZ") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[370,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(as.double(raw_genesE2[994,])) ~ 1 + . , data = covar_Exprsf2[c(1:2,4:5,7:11)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp[370,])) ~ 1 + . + plate, data = covar_Protf[c(1:2,4:5,7)])), BMI =covar_Exprsf2$BMI_cat)
pdf(file = "FAM63A_METRNL_BMI.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe2, y = MolPhe1, color  = BMI)) + geom_point(alpha=0.5,shape=20) +
  scale_color_manual(name = "BMI", values=setNames(c("darkslateblue", "darkslategray", "goldenrod4", "maroon4"),c("Underweight", "Normal", "Overweight", "Obese"))) +
  stat_smooth(data=subset(data, BMI == "Normal"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslategray") +
  #stat_smooth(data=subset(data, BMI == "Overweight"), method = "lm", se = FALSE, fullrange = TRUE, col = "goldenrod4") +
  stat_smooth(data=subset(data, BMI == "Obese"), method = "lm", se = FALSE, fullrange = TRUE, col = "maroon4") +
  labs(title="FAM63A gene expression VS METRNL levels", x ="METRNL", y = "FAM63A") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[221,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(as.double(raw_genesE2[1974,])) ~ 1 + . , data = covar_Exprsf2[c(1:2,4:5,7:11)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp[221,])) ~ 1 + . + plate, data = covar_Protf[c(1:2,4:5,7)])), BMI =covar_Exprsf2$BMI_cat)
pdf(file = "RAB1A_IGFBP2_BMI.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe2, y = MolPhe1, color  = BMI)) + geom_point(alpha=0.5,shape=20) +
  scale_color_manual(name = "BMI", values=setNames(c("darkslateblue", "darkslategray", "goldenrod4", "maroon4"),c("Underweight", "Normal", "Overweight", "Obese"))) +
  stat_smooth(data=subset(data, BMI == "Normal"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslategray") +
  #stat_smooth(data=subset(data, BMI == "Overweight"), method = "lm", se = FALSE, fullrange = TRUE, col = "goldenrod4") +
  stat_smooth(data=subset(data, BMI == "Obese"), method = "lm", se = FALSE, fullrange = TRUE, col = "maroon4") +
  labs(title="RAB1A gene expression VS IGFBP2 levels", x ="IGFBP2", y = "RAB1A") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[172,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(t(raw_met_targE2[1,])) ~ 1 + . , data = covar_Targf2[c(1:2,4:5,7:9)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp[172,])) ~ 1 + . + plate, data = covar_Protf[c(1:2,4:5,7)])), BMI =covar_Exprsf2$BMI_cat)
pdf(file = "GDF15_Carnitine_BMI.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe1, y = MolPhe2, color  = BMI)) + geom_point(alpha=0.5,shape=20) +
  scale_color_manual(name = "BMI", values=setNames(c("darkslateblue", "darkslategray", "goldenrod4", "maroon4"),c("Underweight", "Normal", "Overweight", "Obese"))) +
  stat_smooth(data=subset(data, BMI == "Normal"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslategray") +
  #stat_smooth(data=subset(data, BMI == "Overweight"), method = "lm", se = FALSE, fullrange = TRUE, col = "goldenrod4") +
  stat_smooth(data=subset(data, BMI == "Obese"), method = "lm", se = FALSE, fullrange = TRUE, col = "maroon4") +
  labs(title="GDF15 levels VS L-Carnitine levels", x ="L-Carnitine", y = "GDF15") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf_HbA1cSub[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[79,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(as.double(raw_genesE2_HbA1cSub[994,])) ~ 1 + . , data = covar_Exprsf2_HbA1cSub[c(1:3,5,7:11)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp_HbA1cSub[79,])) ~ 1 + . + plate, data = covar_Protf_HbA1cSub[c(1:3,5,7)])), HbA1c =covar_Exprsf2_HbA1cSub$HbA1c_cat)
pdf(file = "FAM63A_IL6_HbA1c.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe2, y = MolPhe1, color  = HbA1c)) + geom_point(alpha=0.4,shape=20) + scale_color_manual(values=c("maroon4", "darkslategray")) +
  stat_smooth(data=subset(data, HbA1c == "High"), method = "lm", se = FALSE, fullrange = TRUE, col = "maroon4") +
  stat_smooth(data=subset(data, HbA1c == "Low"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslategray") +
  labs(title="FAM63A gene expression VS IL6 levels", x ="IL6", y = "FAM63A") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf_HbA1cSub[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[221,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(t(raw_met_targE2_HbA1cSub[1,])) ~ 1 + . , data = covar_Targf2_AgeSub[c(1:3,5,7:9)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp_HbA1cSub[221,])) ~ 1 + . + plate, data = covar_Protf_HbA1cSub[c(1:3,5,7)])), HbA1c =covar_Exprsf2_HbA1cSub$HbA1c_cat)
pdf(file = "IGFBP2_Carnitine_HbA1c.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe2, y = MolPhe1, color  = HbA1c)) + geom_point(alpha=0.4,shape=20) + scale_color_manual(values=c("maroon4", "darkslategray")) +
  stat_smooth(data=subset(data, HbA1c == "High"), method = "lm", se = FALSE, fullrange = TRUE, col = "maroon4") +
  stat_smooth(data=subset(data, HbA1c == "Low"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslategray") +
  labs(title="IGFBP2 levels VS L-Carnitine levels", x ="L-Carnitine", y = "IGFBP2") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf_HbA1cSub[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[79,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(as.double(raw_genesE2_HbA1cSub[2571,])) ~ 1 + . , data = covar_Exprsf2_HbA1cSub[c(1:3,5,7:11)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp_HbA1cSub[79,])) ~ 1 + . + plate, data = covar_Protf_HbA1cSub[c(1:3,5,7)])), HbA1c =covar_Exprsf2_HbA1cSub$HbA1c_cat)
pdf(file = "IGFBP2_IL6_HbA1c.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe2, y = MolPhe1, color  = HbA1c)) + geom_point(alpha=0.4,shape=20) + scale_color_manual(values=c("maroon4", "darkslategray")) +
  stat_smooth(data=subset(data, HbA1c == "High"), method = "lm", se = FALSE, fullrange = TRUE, col = "maroon4") +
  stat_smooth(data=subset(data, HbA1c == "Low"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslategray") +
  labs(title="IGFBP2 gene expression VS IL6 levels", x ="IL6", y = "IGFBP2") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf_AgeSub[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[172,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(t(raw_met_targE2_AgeSub[1,])) ~ 1 + . , data = covar_Targf2_AgeSub[c(2:5,7:9)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp_AgeSub[172,])) ~ 1 + . + plate, data = covar_Protf_AgeSub[c(2:3,5,7)])), Age_group =covar_Exprsf2_AgeSub$Age_group)
pdf(file = "GDF15_Carnitine_Age.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe1, y = MolPhe2, color  = Age_group)) + geom_point(alpha=0.4,shape=20) + scale_color_manual(values=c("goldenrod1", "darkslateblue")) +
  stat_smooth(data=subset(data, Age_group == "Older"), method = "lm", se = FALSE, fullrange = TRUE, col = "goldenrod1") +
  stat_smooth(data=subset(data, Age_group == "Younger"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslateblue") +
  labs(title="GDF15 levels VS L-Carnitine levels", x ="L-Carnitine", y = "GDF15") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()

plate <- as.factor(as.matrix(covar_Protf_AgeSub[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[370,4],4],3], sep='')]))
data <- data.frame(MolPhe1 = resid(lm(as.matrix(as.double(raw_genesE2_AgeSub[1974,])) ~ 1 + . , data = covar_Exprsf2_AgeSub[c(1:3,5,7:11)])), MolPhe2 = resid(lm(as.matrix(as.double(raw_proteins_Imp_AgeSub[370,])) ~ 1 + . + plate, data = covar_Protf_AgeSub[c(2:3,5,7)])), Age_group =covar_Exprsf2_AgeSub$Age_group)
pdf(file = "RAB1A_METRNL_Age.pdf",  width = 8.3, height = 5.8)
ggplot(data, aes(x = MolPhe2, y = MolPhe1, color  = Age_group)) + geom_point(alpha=0.4,shape=20) + scale_color_manual(values=c("goldenrod1", "darkslateblue")) +
  stat_smooth(data=subset(data, Age_group == "Older"), method = "lm", se = FALSE, fullrange = TRUE, col = "goldenrod1") +
  stat_smooth(data=subset(data, Age_group == "Younger"), method = "lm", se = FALSE, fullrange = TRUE, col = "darkslateblue") +
  labs(title="RAB1A gene expression VS METRNL levels", x ="METRNL", y = "RAB1A") + theme(aspect.ratio=1,title = element_text( size = 16, face = "italic" ),axis.text = element_text(size = 10),axis.title = element_text( size = 16, face = "bold" ),strip.text = element_text(size = 12, face = "bold" ))
dev.off()