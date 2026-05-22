# THIS FILE CONTAINS THE R CODE USED TO IDENTIFY CONTEXT-DEPENDENT EFFECTS ON MOLECULAR PHENOTYPES INCLUDING:
# AGE*SEX, AGE*BMI, AGE*HbA1c, SEX*HbA1c, SEX*BMI, BMI*HbA1c
# FOR GENES, PROTEINS AND METABOLITES (TARGETED & UNTARGETED DATASET)
# USING INTERACTION MODELS
# AND APPLYING MULTIPLE CORRECTION METHODS

#-------------------------------------------------------------------	LIBRARIES/PACKAGES USED	-------------------------------------------------------------------------
library(tibble)
library(dplyr)
library(qvalue)
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#	CONTEXT-DEPENDENT EFFECTS ON MOLECULAR PHENOTYPES: INTERACTION MODELS - FUNCTION

lmCovInter <- function(molPhe, annot, cov1, cov2, otherCov, cov1Name, cov2Name, PROT=FALSE, protein_info=NULL, covariates_prot=NULL){

	results <- matrix(0,ncol=ncol(annot)+6,nrow=nrow(molPhe))
	colnames(results) <- c(colnames(annot), cov1Name,cov2Name,paste(cov1Name,cov2Name,sep="*"),paste(cov1Name,"pvalue",sep="_"),paste(cov2Name,"pvalue",sep="_"),paste(paste(cov1Name,cov2Name,sep="*"),"pvalue",sep="_"))
	panel <- rbind(c('Olink_CARDIOMETABOLIC','C-Met','CAM'), c('Olink_CARDIOVASCULAR_II','CVD_II','CVD2'),
					c('Olink_CARDIOVASCULAR_III','CVD_III','CVD3'), c('Olink_DEVELOPMENT','Dev','DEV'),
					c('Olink_METABOLISM','Met','MET'))

	if(PROT==TRUE){
		for(i in 1:nrow(molPhe)){
			exprsS <- as.matrix(as.double(molPhe[i,]))
			prot <- annot[i,4]
			plateCol <- paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% prot,4],3], sep='')
			plate <- as.factor(as.matrix(covariates_prot[,plateCol]))
			lmInt <- lm(exprsS ~ cov1 * cov2 + cov1 + cov2 + plate + ., data = otherCov)
			coeff <- summary(lmInt)$coefficients
			results[i,] <- c(as.character(annot[i,]),coeff[c(2:3,nrow(coeff)),1],coeff[c(2:3,nrow(coeff)),4])
		}
	}	
	else{
		for(i in 1:nrow(molPhe)){
			exprsS <- as.matrix(as.double(molPhe[i,]))
			lmInt <- lm(exprsS ~ cov1 * cov2 + cov1 + cov2 + ., data = otherCov)
			coeff <- summary(lmInt)$coefficients
			results[i,] <- c(as.character(annot[i,]),coeff[c(2:3,nrow(coeff)),1],coeff[c(2:3,nrow(coeff)),4])
		}
	}
	return(results)
}


#		CONTEXT-DEPENDENT EFFECTS ON MOLECULAR PHENOTYPES: INTERACTION MODELS - GET RESULTS
##	Scale Covariates
covar_Exprsf_S <- covar_Exprsf
covar_Exprsf_S$Age <- scale(covar_Exprsf$Age, scale = FALSE)
covar_Exprsf_S$BMI <- scale(covar_Exprsf$BMI, scale = FALSE)
covar_Exprsf_S$GC_mean <- scale(covar_Exprsf$GC_mean, scale=FALSE)
covar_Exprsf_S$InsertSizeMode <- scale(covar_Exprsf$InsertSizeMode, scale=FALSE)
covar_Exprsf_S$Rin <- scale(covar_Exprsf$Rin, scale=FALSE)
covar_Exprsf_S$HbA1c <- scale(covar_Exprsf$HbA1c , scale=FALSE)

covar_Protf_S <- covar_Protf
covar_Protf_S$Age <- scale(covar_Protf$Age, scale = FALSE)
covar_Protf_S$BMI <- scale(covar_Protf$BMI, scale = FALSE)
covar_Protf_S$HbA1c <- scale(covar_Protf$HbA1c , scale=FALSE)

covar_Targf_S <- covar_Targf
covar_Targf_S$Age <- scale(covar_Targf$Age, scale = FALSE)
covar_Targf_S$BMI <- scale(covar_Targf$BMI, scale = FALSE)
covar_Targf_S$HbA1c <- scale(covar_Targf$HbA1c , scale=FALSE)

covar_Untargf_S <- covar_Untargf
covar_Untargf_S$Age <- scale(covar_Untargf$Age, scale = FALSE)
covar_Untargf_S$BMI <- scale(covar_Untargf$BMI, scale = FALSE)
covar_Untargf_S$HbA1c <- scale(covar_Untargf$HbA1c , scale=FALSE)

#*****************************************************************************************		GENES		**************************************************************************************************************
res_Genes_AgeSex <- lmCovInter(
										raw_genesE,
										annot=anno_genes,
										cov1=covar_Exprsf_S$Age,
										cov2=covar_Exprsf_S$Sex,
										otherCov = covar_Exprsf_S[c(3:5,7:11)],
										cov1Name="Age",
										cov2Name="Sex")
res_Genes_AgeBMI <- lmCovInter(
										raw_genesE,
										annot=anno_genes,
										cov1=covar_Exprsf_S$Age,
										cov2=covar_Exprsf_S$BMI,
										otherCov = covar_Exprsf_S[c(2,4:5,7:11)],
										cov1Name="Age",
										cov2Name="BMI")
res_Genes_AgeHbA1c <- lmCovInter(
											raw_genesE,
											annot=anno_genes,
											cov1=covar_Exprsf_S$Age,
											cov2=covar_Exprsf_S$HbA1c,
											otherCov = covar_Exprsf_S[c(2:3,5,7:11)],
											cov1Name="Age",
											cov2Name="HbA1c")
res_Genes_SexHbA1c <- lmCovInter(
											raw_genesE,
											annot=anno_genes,
											cov1=covar_Exprsf_S$Sex,
											cov2=covar_Exprsf_S$HbA1c,
											otherCov = covar_Exprsf_S[c(1,3,5,7:11)],
											cov1Name="Sex",
											cov2Name="HbA1c")
res_Genes_SexBMI <- lmCovInter(
										raw_genesE,
										annot=anno_genes,
										cov1=covar_Exprsf_S$Sex,
										cov2=covar_Exprsf_S$BMI,
										otherCov = covar_Exprsf_S[c(1,4:5,7:11)],
										cov1Name="Sex",
										cov2Name="BMI")
res_Genes_BMIHbA1c <- lmCovInter(
											raw_genesE,
											annot=anno_genes,
											cov1=covar_Exprsf_S$BMI,
											cov2=covar_Exprsf_S$HbA1c,
											otherCov = covar_Exprsf_S[c(1:2,5,7:11)],
											cov1Name="BMI",
											cov2Name="HbA1c")
res_Genes_AgeSex <- as.data.frame(res_Genes_AgeSex)
res_Genes_AgeBMI <- as.data.frame(res_Genes_AgeBMI)
res_Genes_AgeHbA1c <- as.data.frame(res_Genes_AgeHbA1c)
res_Genes_SexHbA1c <- as.data.frame(res_Genes_SexHbA1c)
res_Genes_SexBMI <- as.data.frame(res_Genes_SexBMI)
res_Genes_BMIHbA1c <- as.data.frame(res_Genes_BMIHbA1c)
#****************************************************************************************			PROTEINS			*****************************************************************************************************
res_Proteins_AgeSex <- lmCovInter(
											raw_proteins_Imp,
											annot=anno_proteins,
											cov1=covar_Protf_S$Age,
											cov2=covar_Protf_S$Sex,
											otherCov = covar_Protf_S[c(3:5,7)],
											cov1Name="Age",
											cov2Name="Sex",
											PROT=TRUE,
											protein_info=protein_info,
											covariates_prot=covar_Protf_S[8:12])
res_Proteins_AgeBMI <- lmCovInter(
											raw_proteins_Imp,
											annot=anno_proteins,
											cov1=covar_Protf_S$Age,
											cov2=covar_Protf_S$BMI,
											otherCov = covar_Protf_S[c(2,4:5,7)],
											cov1Name="Age",
											cov2Name="BMI",
											PROT=TRUE,
											protein_info=protein_info,
											covariates_prot=covar_Protf_S[8:12])
res_Proteins_AgeHbA1c <- lmCovInter(
												raw_proteins_Imp,
												annot=anno_proteins,
												cov1=covar_Protf_S$Age,
												cov2=covar_Protf_S$HbA1c,
												otherCov = covar_Protf_S[c(2:5,7)],
												cov1Name="Age",
												cov2Name="HbA1c",
												PROT=TRUE,
												protein_info=protein_info,
											covariates_prot=covar_Protf_S[8:12])
res_Proteins_SexHbA1c <- lmCovInter(
											raw_proteins_Imp,
											annot=anno_proteins,
											cov1=covar_Protf_S$Sex,
											cov2=covar_Protf_S$HbA1c,
											otherCov = covar_Protf_S[c(1,3,5,7)],
											cov1Name="Sex",
											cov2Name="HbA1c",
											PROT=TRUE,
											protein_info=protein_info,
											covariates_prot=covar_Protf_S[8:12])
res_Proteins_SexBMI <- lmCovInter(
											raw_proteins_Imp,
											annot=anno_proteins,
											cov1=covar_Protf_S$Sex,
											cov2=covar_Protf_S$BMI,
											otherCov = covar_Protf_S[c(1,4:5,7)],
											cov1Name="Sex",
											cov2Name="BMI",
											PROT=TRUE,
											protein_info=protein_info,
											covariates_prot=covar_Protf_S[8:12])
res_Proteins_BMIHbA1c <- lmCovInter(
											raw_proteins_Imp,
											annot=anno_proteins,
											cov1=covar_Protf_S$BMI,
											cov2=covar_Protf_S$HbA1c,
											otherCov = covar_Protf_S[c(1:2,5,7)],
											cov1Name="BMI",
											cov2Name="HbA1c",
											PROT=TRUE,
											protein_info=protein_info,
											covariates_prot=covar_Protf_S[8:12])													
res_Proteins_AgeSex <- as.data.frame(res_Proteins_AgeSex)
res_Proteins_AgeBMI <- as.data.frame(res_Proteins_AgeBMI)
res_Proteins_AgeHbA1c <- as.data.frame(res_Proteins_AgeHbA1c)
res_Proteins_SexHbA1c <- as.data.frame(res_Proteins_SexHbA1c)
res_Proteins_SexBMI <- as.data.frame(res_Proteins_SexBMI)
res_Proteins_BMIHbA1c <- as.data.frame(res_Proteins_BMIHbA1c)
#********************************************************************************			TARGETED METABOLITES			********************************************************************************************
res_Targ_AgeSex <- lmCovInter(
										raw_met_targE,
										annot=anno_met_targ,
										cov1=covar_Targf_S$Age,
										cov2=covar_Targf_S$Sex,
										otherCov = covar_Targf_S[c(3:5,7:9)],
										cov1Name="Age",
										cov2Name="Sex")
res_Targ_AgeBMI <- lmCovInter(
										raw_met_targE,
										annot=anno_met_targ,
										cov1=covar_Targf_S$Age,
										cov2=covar_Targf_S$BMI,
										otherCov = covar_Targf_S[c(2,4:5,7:9)],
										cov1Name="Age",
										cov2Name="BMI")
res_Targ_AgeHbA1c <- lmCovInter(
											raw_met_targE,
											annot=anno_met_targ,
											cov1=covar_Targf_S$Age,
											cov2=covar_Targf_S$HbA1c,
											otherCov = covar_Targf_S[c(2:3,5,7:9)],
											cov1Name="Age",
											cov2Name="HbA1c")	
res_Targ_SexHbA1c <- lmCovInter(
											raw_met_targE,
											annot=anno_met_targ,
											cov1=covar_Targf_S$Sex,
											cov2=covar_Targf_S$HbA1c,
											otherCov = covar_Targf_S[c(1,3,5,7:9)],
											cov1Name="Sex",
											cov2Name="HbA1c")
res_Targ_SexBMI <- lmCovInter(
										raw_met_targE,
										annot=anno_met_targ,
										cov1=covar_Targf_S$Sex,
										cov2=covar_Targf_S$BMI,
										otherCov = covar_Targf_S[c(1,4:5,7:9)],
										cov1Name="Sex",
										cov2Name="BMI")
res_Targ_BMIHbA1c <- lmCovInter(
										raw_met_targE,
										annot=anno_met_targ,
										cov1=covar_Targf_S$BMI,
										cov2=covar_Targf_S$HbA1c,
										otherCov = covar_Targf_S[c(1:2,5,7:9)],
										cov1Name="BMI",
										cov2Name="HbA1c")											
res_Targ_AgeSex <- as.data.frame(res_Targ_AgeSex)
res_Targ_AgeBMI <- as.data.frame(res_Targ_AgeBMI)
res_Targ_AgeHbA1c <- as.data.frame(res_Targ_AgeHbA1c)
res_Targ_SexHbA1c <- as.data.frame(res_Targ_SexHbA1c)
res_Targ_SexBMI <- as.data.frame(res_Targ_SexBMI)
res_Targ_BMIHbA1c <- as.data.frame(res_Targ_BMIHbA1c)
#********************************************************************************			UNTARGETED METABOLITES			********************************************************************************************
res_Untarg_AgeSex <- lmCovInter(
										raw_met_untarg_Imp,
										annot=anno_met_untarg,
										cov1=covar_Untargf_S$Age,
										cov2=covar_Untargf_S$Sex,
										otherCov = covar_Untargf_S[c(3:5,7:8)],
										cov1Name="Age",
										cov2Name="Sex")
res_Untarg_AgeBMI <- lmCovInter(
										raw_met_untarg_Imp,
										annot=anno_met_untarg,
										cov1=covar_Untargf_S$Age,
										cov2=covar_Untargf_S$BMI,
										otherCov = covar_Untargf_S[c(2,4:5,7:8)],
										cov1Name="Age",
										cov2Name="BMI")
res_Untarg_AgeHbA1c <- lmCovInter(
											raw_met_untarg_Imp,
											annot=anno_met_untarg,
											cov1=covar_Untargf_S$Age,
											cov2=covar_Untargf_S$HbA1c,
											otherCov = covar_Untargf_S[c(2:3,5,7:8)],
											cov1Name="Age",
											cov2Name="HbA1c")
res_Untarg_SexHbA1c <- lmCovInter(
											raw_met_untarg_Imp,
											annot=anno_met_untarg,
											cov1=covar_Untargf_S$Sex,
											cov2=covar_Untargf_S$HbA1c,
											otherCov = covar_Untargf_S[c(1,3,5,7:8)],
											cov1Name="Sex",
											cov2Name="HbA1c")
res_Untarg_SexBMI <- lmCovInter(
										raw_met_untarg_Imp,
										annot=anno_met_untarg,
										cov1=covar_Untargf_S$Sex,
										cov2=covar_Untargf_S$BMI,
										otherCov = covar_Untargf_S[c(1,4:5,7:8)],
										cov1Name="Sex",
										cov2Name="BMI")
res_Untarg_BMIHbA1c <- lmCovInter(
											raw_met_untarg_Imp,
											annot=anno_met_untarg,
											cov1=covar_Untargf_S$BMI,
											cov2=covar_Untargf_S$HbA1c,
											otherCov = covar_Untargf_S[c(1:2,5,7:8)],
											cov1Name="BMI",
											cov2Name="HbA1c")
res_Untarg_AgeSex <- as.data.frame(res_Untarg_AgeSex)
res_Untarg_AgeBMI <- as.data.frame(res_Untarg_AgeBMI)
res_Untarg_AgeHbA1c <- as.data.frame(res_Untarg_AgeHbA1c)
res_Untarg_SexHbA1c <- as.data.frame(res_Untarg_SexHbA1c)
res_Untarg_SexBMI <- as.data.frame(res_Untarg_SexBMI)
res_Untarg_BMIHbA1c <- as.data.frame(res_Untarg_BMIHbA1c)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#	MULTIPLE CORRECTION
res_Genes_AgeSexAdj <- res_Genes_AgeSex
res_Genes_AgeSexAdj$Age_Q <- qvalue(as.numeric(res_Genes_AgeSexAdj[,"Age_pvalue"]))$qvalues
res_Genes_AgeSexAdj$Sex_Q <- qvalue(as.numeric(res_Genes_AgeSexAdj[,"Sex_pvalue"]))$qvalues
res_Genes_AgeSexAdj$"Age*Sex_Q" <- qvalue(as.numeric(res_Genes_AgeSexAdj[,"Age*Sex_pvalue"]))$qvalues
res_Genes_AgeBMIAdj <- res_Genes_AgeBMI
res_Genes_AgeBMIAdj$Age_Q <- qvalue(as.numeric(res_Genes_AgeBMIAdj[,"Age_pvalue"]))$qvalues
res_Genes_AgeBMIAdj$BMI_Q <- qvalue(as.numeric(res_Genes_AgeBMIAdj[,"BMI_pvalue"]))$qvalues
res_Genes_AgeBMIAdj$"Age*BMI_Q" <- qvalue(as.numeric(res_Genes_AgeBMIAdj[,"Age*BMI_pvalue"]))$qvalues
res_Genes_AgeHbA1cAdj <- res_Genes_AgeHbA1c
res_Genes_AgeHbA1cAdj$Age_Q <- qvalue(as.numeric(res_Genes_AgeHbA1cAdj[,"Age_pvalue"]))$qvalues
res_Genes_AgeHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Genes_AgeHbA1cAdj[,"HbA1c_pvalue"]))$qvalues
res_Genes_AgeHbA1cAdj$"Age*HbA1c_Q" <- qvalue(as.numeric(res_Genes_AgeHbA1cAdj[,"Age*HbA1c_pvalue"]))$qvalues
res_Genes_SexHbA1cAdj <- res_Genes_SexHbA1c
res_Genes_SexHbA1cAdj$Sex_Q <- qvalue(as.numeric(res_Genes_SexHbA1cAdj[,"Sex_pvalue"]))$qvalues
res_Genes_SexHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Genes_SexHbA1cAdj[,"HbA1c_pvalue"]))$qvalues
res_Genes_SexHbA1cAdj$"Sex*HbA1c_Q" <- qvalue(as.numeric(res_Genes_SexHbA1cAdj[,"Sex*HbA1c_pvalue"]))$qvalues
res_Genes_SexBMIAdj <- res_Genes_SexBMI
res_Genes_SexBMIAdj$Sex_Q <- qvalue(as.numeric(res_Genes_SexBMIAdj[,"Sex_pvalue"]))$qvalues
res_Genes_SexBMIAdj$BMI_Q <- qvalue(as.numeric(res_Genes_SexBMIAdj[,"BMI_pvalue"]))$qvalues
res_Genes_SexBMIAdj$"Sex*BMI_Q" <- qvalue(as.numeric(res_Genes_SexBMIAdj[,"Sex*BMI_pvalue"]))$qvalues
res_Genes_BMIHbA1cAdj <- res_Genes_BMIHbA1c
res_Genes_BMIHbA1cAdj$BMI_Q <- qvalue(as.numeric(res_Genes_BMIHbA1cAdj[,"BMI_pvalue"]))$qvalues
res_Genes_BMIHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Genes_BMIHbA1cAdj[,"HbA1c_pvalue"]))$qvalues
res_Genes_BMIHbA1cAdj$"BMI*HbA1c_Q" <- qvalue(as.numeric(res_Genes_BMIHbA1cAdj[,"BMI*HbA1c_pvalue"]))$qvalues

res_Proteins_AgeSexAdj <- res_Proteins_AgeSex
res_Proteins_AgeSexAdj$Age_Q <- qvalue(as.numeric(res_Proteins_AgeSexAdj[,"Age_pvalue"]))$qvalues
res_Proteins_AgeSexAdj$Sex_Q <- qvalue(as.numeric(res_Proteins_AgeSexAdj[,"Sex_pvalue"]))$qvalues
res_Proteins_AgeSexAdj$"Age*Sex_Q" <- qvalue(as.numeric(res_Proteins_AgeSexAdj[,"Age*Sex_pvalue"]))$qvalues
res_Proteins_AgeBMIAdj <- res_Proteins_AgeBMI
res_Proteins_AgeBMIAdj$Age_Q <- qvalue(as.numeric(res_Proteins_AgeBMIAdj[,"Age_pvalue"]))$qvalues
res_Proteins_AgeBMIAdj$BMI_Q <- qvalue(as.numeric(res_Proteins_AgeBMIAdj[,"BMI_pvalue"]))$qvalues
res_Proteins_AgeBMIAdj$"Age*BMI_Q" <- qvalue(as.numeric(res_Proteins_AgeBMIAdj[,"Age*BMI_pvalue"]))$qvalues
res_Proteins_AgeHbA1cAdj <- res_Proteins_AgeHbA1c
res_Proteins_AgeHbA1cAdj$Age_Q <- qvalue(as.numeric(res_Proteins_AgeHbA1cAdj[,"Age_pvalue"]))$qvalues
res_Proteins_AgeHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Proteins_AgeHbA1cAdj[,"HbA1c_pvalue"]))$qvalues
res_Proteins_AgeHbA1cAdj$"Age*HbA1c_Q" <- qvalue(as.numeric(res_Proteins_AgeHbA1cAdj[,"Age*HbA1c_pvalue"]))$qvalues
res_Proteins_SexHbA1cAdj <- res_Proteins_SexHbA1c
res_Proteins_SexHbA1cAdj$Sex_Q <- qvalue(as.numeric(res_Proteins_SexHbA1cAdj[,"Sex_pvalue"]))$qvalues
res_Proteins_SexHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Proteins_SexHbA1cAdj[,"HbA1c_pvalue"]))$qvalues
res_Proteins_SexHbA1cAdj$"Sex*HbA1c_Q" <- qvalue(as.numeric(res_Proteins_SexHbA1cAdj[,"Sex*HbA1c_pvalue"]))$qvalues
res_Proteins_SexBMIAdj <- res_Proteins_SexBMI
res_Proteins_SexBMIAdj$Sex_Q <- qvalue(as.numeric(res_Proteins_SexBMIAdj[,"Sex_pvalue"]))$qvalues
res_Proteins_SexBMIAdj$BMI_Q <- qvalue(as.numeric(res_Proteins_SexBMIAdj[,"BMI_pvalue"]))$qvalues
res_Proteins_SexBMIAdj$"Sex*BMI_Q" <- qvalue(as.numeric(res_Proteins_SexBMIAdj[,"Sex*BMI_pvalue"]))$qvalues
res_Proteins_BMIHbA1cAdj <- res_Proteins_BMIHbA1c
res_Proteins_BMIHbA1cAdj$BMI_Q <- qvalue(as.numeric(res_Proteins_BMIHbA1cAdj[,"BMI_pvalue"]), lambda = seq(0, max(as.numeric(res_Proteins_BMIHbA1cAdj$BMI_pvalue)), 0.05))$qvalues
res_Proteins_BMIHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Proteins_BMIHbA1cAdj[,"HbA1c_pvalue"]))$qvalues
res_Proteins_BMIHbA1cAdj$"BMI*HbA1c_Q" <- qvalue(as.numeric(res_Proteins_BMIHbA1cAdj[,"BMI*HbA1c_pvalue"]))$qvalues

res_Targ_AgeSexAdj <- res_Targ_AgeSex
res_Targ_AgeSexAdj$Age_Q <- qvalue(as.numeric(res_Targ_AgeSexAdj[,"Age_pvalue"]))$qvalues
res_Targ_AgeSexAdj$Sex_Q <- qvalue(as.numeric(res_Targ_AgeSexAdj[,"Sex_pvalue"]), lambda=c(0.01,0.05,0.1,0.5,0.9))$qvalues
res_Targ_AgeSexAdj$"Age*Sex_Q" <- qvalue(as.numeric(res_Targ_AgeSexAdj[,"Age*Sex_pvalue"]))$qvalues
res_Targ_AgeBMIAdj <- res_Targ_AgeBMI
res_Targ_AgeBMIAdj$Age_Q <- qvalue(as.numeric(res_Targ_AgeBMIAdj[,"Age_pvalue"]), lambda=c(0.01,0.05,0.1,0.5,0.9))$qvalues
res_Targ_AgeBMIAdj$BMI_Q <- qvalue(as.numeric(res_Targ_AgeBMIAdj[,"BMI_pvalue"]))$qvalues
res_Targ_AgeBMIAdj$"Age*BMI_Q" <- qvalue(as.numeric(res_Targ_AgeBMIAdj[,"Age*BMI_pvalue"]))$qvalues
res_Targ_AgeHbA1cAdj <- res_Targ_AgeHbA1c
res_Targ_AgeHbA1cAdj$Age_Q <- qvalue(as.numeric(res_Targ_AgeHbA1cAdj[,"Age_pvalue"]), lambda = seq(0, max(as.numeric(res_Targ_AgeHbA1cAdj[,"Age_pvalue"])), 0.05))$qvalues
res_Targ_AgeHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Targ_AgeHbA1cAdj[,"HbA1c_pvalue"]), lambda=c(0.01,0.05,0.1,0.5,0.9))$qvalues
res_Targ_AgeHbA1cAdj$"Age*HbA1c_Q" <- qvalue(as.numeric(res_Targ_AgeHbA1cAdj[,"Age*HbA1c_pvalue"]))$qvalues
res_Targ_SexHbA1cAdj <- res_Targ_SexHbA1c
res_Targ_SexHbA1cAdj$Sex_Q <- qvalue(as.numeric(res_Targ_SexHbA1cAdj[,"Sex_pvalue"]))$qvalues
res_Targ_SexHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Targ_SexHbA1cAdj[,"HbA1c_pvalue"]))$qvalues
res_Targ_SexHbA1cAdj$"Sex*HbA1c_Q" <- qvalue(as.numeric(res_Targ_SexHbA1cAdj[,"Sex*HbA1c_pvalue"]))$qvalues
res_Targ_SexBMIAdj <- res_Targ_SexBMI
res_Targ_SexBMIAdj$Sex_Q <- qvalue(as.numeric(res_Targ_SexBMIAdj[,"Sex_pvalue"]))$qvalues
res_Targ_SexBMIAdj$BMI_Q <- qvalue(as.numeric(res_Targ_SexBMIAdj[,"BMI_pvalue"]))$qvalues
res_Targ_SexBMIAdj$"Sex*BMI_Q" <- qvalue(as.numeric(res_Targ_SexBMIAdj[,"Sex*BMI_pvalue"]))$qvalues
res_Targ_BMIHbA1cAdj <- res_Targ_BMIHbA1c
res_Targ_BMIHbA1cAdj$BMI_Q <- qvalue(as.numeric(res_Targ_BMIHbA1cAdj[,"BMI_pvalue"]))$qvalues
res_Targ_BMIHbA1cAdj$HbA1c_Q <-  qvalue(as.numeric(res_Targ_BMIHbA1cAdj[,"HbA1c_pvalue"]), lambda=c(0.01,0.05,0.1,0.5,0.9))$qvalues
res_Targ_BMIHbA1cAdj$"BMI*HbA1c_Q" <- qvalue(as.numeric(res_Targ_BMIHbA1cAdj[,"BMI*HbA1c_pvalue"]))$qvalues

res_Untarg_AgeSexAdj <- res_Untarg_AgeSex
res_Untarg_AgeSexAdj$Age_Q <- qvalue(as.numeric(res_Untarg_AgeSexAdj[,"Age_pvalue"]))$qvalues
res_Untarg_AgeSexAdj$Sex_Q <- qvalue(as.numeric(res_Untarg_AgeSexAdj[,"Sex_pvalue"]))$qvalues
res_Untarg_AgeSexAdj$"Age*Sex_Q" <- qvalue(as.numeric(res_Untarg_AgeSexAdj[,"Age*Sex_pvalue"]))$qvalues
res_Untarg_AgeBMIAdj <- res_Untarg_AgeBMI
res_Untarg_AgeBMIAdj$Age_Q <- qvalue(as.numeric(res_Untarg_AgeBMIAdj[,"Age_pvalue"]))$qvalues
res_Untarg_AgeBMIAdj$BMI_Q <- qvalue(as.numeric(res_Untarg_AgeBMIAdj[,"BMI_pvalue"]), lambda=c(0.01,0.05,0.1,0.5,0.9))$qvalues
res_Untarg_AgeBMIAdj$"Age*BMI_Q" <- qvalue(as.numeric(res_Untarg_AgeBMIAdj[,"Age*BMI_pvalue"]))$qvalues
res_Untarg_AgeHbA1cAdj <- res_Untarg_AgeHbA1c
res_Untarg_AgeHbA1cAdj$Age_Q <- qvalue(as.numeric(res_Untarg_AgeHbA1cAdj[,"Age_pvalue"]))$qvalues
res_Untarg_AgeHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Untarg_AgeHbA1cAdj[,"HbA1c_pvalue"]))$qvalues
res_Untarg_AgeHbA1cAdj$"Age*HbA1c_Q" <- qvalue(as.numeric(res_Untarg_AgeHbA1cAdj[,"Age*HbA1c_pvalue"]))$qvalues
res_Untarg_SexHbA1cAdj <- res_Untarg_SexHbA1c
res_Untarg_SexHbA1cAdj$Sex_Q <- qvalue(as.numeric(res_Untarg_SexHbA1cAdj[,"Sex_pvalue"]), lambda = seq(0, max(as.numeric(res_Untarg_SexHbA1cAdj[,"Sex_pvalue"])), 0.05))$qvalues
res_Untarg_SexHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Untarg_SexHbA1cAdj[,"HbA1c_pvalue"]))$qvalues
res_Untarg_SexHbA1cAdj$"Sex*HbA1c_Q" <- qvalue(as.numeric(res_Untarg_SexHbA1cAdj[,"Sex*HbA1c_pvalue"]))$qvalues
res_Untarg_SexBMIAdj <- res_Untarg_SexBMI
res_Untarg_SexBMIAdj$Sex_Q <- qvalue(as.numeric(res_Untarg_SexBMIAdj[,"Sex_pvalue"]))$qvalues
res_Untarg_SexBMIAdj$BMI_Q <- qvalue(as.numeric(res_Untarg_SexBMIAdj[,"BMI_pvalue"]))$qvalues
res_Untarg_SexBMIAdj$"Sex*BMI_Q" <- qvalue(as.numeric(res_Untarg_SexBMIAdj[,"Sex*BMI_pvalue"]))$qvalues
res_Untarg_BMIHbA1cAdj <- res_Untarg_BMIHbA1c
res_Untarg_BMIHbA1cAdj$BMI_Q <- qvalue(as.numeric(res_Untarg_BMIHbA1cAdj[,"BMI_pvalue"]))$qvalues
res_Untarg_BMIHbA1cAdj$HbA1c_Q <- qvalue(as.numeric(res_Untarg_BMIHbA1cAdj[,"HbA1c_pvalue"]))$qvalues
res_Untarg_BMIHbA1cAdj$"BMI*HbA1c_Q" <- qvalue(as.numeric(res_Untarg_BMIHbA1cAdj[,"BMI*HbA1c_pvalue"]))$qvalues
