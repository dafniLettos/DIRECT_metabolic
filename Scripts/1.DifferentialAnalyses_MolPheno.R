# THIS FILE CONTAINS THE R CODE USED TO IDENTIFY MOLECULAR PHENOTYPES INCLUDING:
# GENES, PROTEINS AND METABOLITES (TARGETED & UNTARGETED DATASET)
# ASSOCIATED WITH: (A).AGE	(B).SEX	(C).BMI	(D).HBA1C
# USING LINEAR REGRESSION MODELS
# AND APPLYING MULTIPLE CORRECTION METHODS
# DONE FOR DIRECT BASELINE: 3,027 INDIVIDUALS

#-------------------------------------------------------------------	LIBRARIES/PACKAGES USED	-------------------------------------------------------------------------
library(tibble)
library(dplyr)
library(qvalue)
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEFINE FUNCTIONS
`%!in%` <- Negate(`%in%`)

normal <- function(v){
	w <- rep(0, length(v))
    ranking <- rank(v, ties.method = "min")
    n <- length(unique(ranking))
    w[1:length(v) %in% ranking] <- qnorm(seq(from = 1/(n + 1), to = n/(n + 1), by = 1/(n + 1)))
    return(w[ranking])
}

#	DATA PRE-PROCESSING
##	Impute Missing Values
raw_molPhenoE <- as.matrix(raw_molPheno[,7:ncol(raw_molPheno)])
names(raw_molPhenoE) <- NULL
raw_molPheno_Imp <- raw_molPhenoE
for(i in 1:nrow(raw_molPheno_Imp)){
	k <- which(is.na(raw_molPheno_Imp[i,]), arr.ind=TRUE)
	for(j in k){
		raw_molPheno_Imp[i,j] <- mean(raw_molPheno_Imp[i,], na.rm = TRUE)}}

##	Rank Normalisation
rownames(raw_molPheno_Imp) <- 1:nrow(raw_molPheno_Imp)
for(i in 1:nrow(raw_molPheno_Imp)){raw_proteins_Imp[i,] <-  normal(as.numeric(raw_molPheno_Imp[i,]))}
##	Test Normalisation
KS_prot <- matrix(0,ncol=3,nrow=nrow(raw_molPheno_Imp))
colnames(KS_prot) <- c("MolPheno", "ks_Value", "ks_Pval")
KS_prot[,1] <- raw_molPheno_Imp$MolPhenoID
for(i in 1:nrow(KS_prot)){
	KS_prot[i,2] <- ks.test(raw_molPheno_Imp[i,],y="pnorm")$statistic
	KS_prot[i,3] <- ks.test(raw_molPheno_Imp[i,],y="pnorm")$p.value}
sum(KS_prot[,3] < 0.05)


#	DIFFERENTIAL ANALYSES: LINEAR MODELS - FUNCTION
lmMolPheno <- function(molPhe, annot, covariates, NoOfCov = 5, PROT=FALSE, protein_info=NULL, covariates_prot=NULL){

	results <- matrix(0,ncol = ncol(annot) + 3*NoOfCov, nrow=nrow(molPhe))
	colnames(results) <- c(colnames(annot), colnames(covariates)[1:NoOfCov], paste(colnames(covariates)[1:NoOfCov],"Tvalue",sep="_"), paste(colnames(covariates)[1:NoOfCov],"Pvalue",sep="_"))
	panel <- rbind(c('Olink_CARDIOMETABOLIC','C-Met','CAM'), c('Olink_CARDIOVASCULAR_II','CVD_II','CVD2'),
					c('Olink_CARDIOVASCULAR_III','CVD_III','CVD3'), c('Olink_DEVELOPMENT','Dev','DEV'),
					c('Olink_METABOLISM','Met','MET'))

	if(PROT==TRUE){
		for(i in 1:nrow(molPhe)){
			molPhe1 <- as.matrix(as.double(molPhe[i,]))
			prot <- annot[i,4]
			plateCol <- paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% prot,4],3], sep='')
			plate <- as.factor(as.matrix(covariates_prot[,plateCol]))
			pdata <- data.frame(covariates, plate)
			lmMolPhe <- lm(molPhe1 ~ 1 + ., data = pdata)
			coeff <- summary(lmMolPhe)$coefficients
			results[i,] <- c(as.character(annot[i,]),coeff[2:(NoOfCov+1),1],coeff[2:(NoOfCov+1),3],coeff[2:(NoOfCov+1),4])
		}
	}	
	else{
		for(i in 1:nrow(molPhe)){
			molPhe1 <- as.matrix(as.double(molPhe[i,]))
			lmMolPhe <- lm(molPhe1 ~ 1 + ., data = covariates)
			coeff <- summary(lmMolPhe)$coefficients
			results[i,] <- c(as.character(annot[i,]),coeff[2:(NoOfCov+1),1],coeff[2:(NoOfCov+1),3],coeff[2:(NoOfCov+1),4])
		}
	}
	return(results)
}


#		DIFFERENTIAL ANALYSES: LINEAR MODELS - GET RESULTS
##	Gene Expression ~ Age + Sex + BMI + HbA1c + WP + Centre + GC Mean Content + Insert Size + Rin Score + Date of Sequencing
res_Genes <- lmPheno(raw_genesE,
									annot = anno_genes,
									covariates = covar_Exprsf[c(1:5,7:11)])
res_Genes <- as.data.frame(res_Genes)
dim(res_Genes)
#[1]	16209 20

##	Protein Levels ~ Age + Sex + BMI + HbA1c + WP + Centre + Protein Plate
res_Proteins <- lmPheno(raw_proteins_Imp,
									annot = anno_proteins,
									covariates = covar_Protf[c(1:5,7)],
									PROT = TRUE,
									protein_info = protein_info,
									covariates_prot = covar_Protf[8:12])
res_Proteins <- as.data.frame(res_Proteins)									
dim(res_Proteins)
#[1]	373 21

##	Metabolite Levels ~ Age + Sex + BMI + HbA1c + WP + Centre + Metabolite Plate + Date of Processing
res_Targ <- lmPheno(raw_met_targE,
								annot = anno_met_targ,
								covariates = covar_Targf[c(1:5,7:9)])
res_Targ <- as.data.frame(res_Targ)
dim(res_Targ)
#[1]	116 16

##	Metabolite Levels ~ Age + Sex + BMI + HbA1c + WP + Centre + Metabolite Plate + Date of Processing
res_Untarg <- lmPheno(raw_met_untarg_Imp,
								annot = anno_met_untarg,
								covariates = covar_Untargf[c(1:5,7:8)])
res_Untarg <- as.data.frame(res_Untarg)
dim(res_Untarg)
#[1]	233 19



#	MULTIPLE CORRECTION
res_GenesAdj <- res_Genes
res_GenesAdj$Age_FDR <- p.adjust(as.numeric(res_GenesAdj$Age_pvalue), method = "BH")
res_GenesAdj$Sex_FDR <- p.adjust(as.numeric(res_GenesAdj$Sex_pvalue), method = "BH")
res_GenesAdj$BMI_FDR <- p.adjust(as.numeric(res_GenesAdj$BMI_pvalue), method = "BH")
res_GenesAdj$HbA1c_FDR <- p.adjust(as.numeric(res_GenesAdj$HbA1c_pvalue), method = "BH")
res_GenesAdj$WP_FDR <- p.adjust(as.numeric(res_GenesAdj$WP_pvalue), method = "BH")
res_GenesAdj$Age_Q <- qvalue(as.numeric(res_GenesAdj$Age_pvalue))$qvalues
res_GenesAdj$Sex_Q <- qvalue(as.numeric(res_GenesAdj$Sex_pvalue))$qvalues
res_GenesAdj$BMI_Q <- qvalue(as.numeric(res_GenesAdj$BMI_pvalue))$qvalues
res_GenesAdj$HbA1c_Q <- qvalue(as.numeric(res_GenesAdj$HbA1c_pvalue))$qvalues
res_GenesAdj$WP_Q <- qvalue(as.numeric(res_GenesAdj$WP_pvalue))$qvalues

res_ProteinsAdj <- res_Proteins
res_ProteinsAdj$Age_FDR <- p.adjust(as.numeric(res_ProteinsAdj$Age_pvalue), method = "BH")
res_ProteinsAdj$Sex_FDR <- p.adjust(as.numeric(res_ProteinsAdj$Sex_pvalue), method = "BH")
res_ProteinsAdj$BMI_FDR <- p.adjust(as.numeric(res_ProteinsAdj$BMI_pvalue), method = "BH")
res_ProteinsAdj$HbA1c_FDR <- p.adjust(as.numeric(res_ProteinsAdj$HbA1c_pvalue), method = "BH")
res_ProteinsAdj$WP_FDR <- p.adjust(as.numeric(res_ProteinsAdj$WP_pvalue), method = "BH")
res_ProteinsAdj$Age_Q <- qvalue(as.numeric(res_ProteinsAdj$Age_pvalue))$qvalues
res_ProteinsAdj$Sex_Q <- qvalue(as.numeric(res_ProteinsAdj$Sex_pvalue))$qvalues
res_ProteinsAdj$BMI_Q <- qvalue(as.numeric(res_ProteinsAdj$BMI_pvalue), lambda = seq(0, max(as.numeric(res_ProteinsAdj$BMI_pvalue)), 0.05))$qvalues
res_ProteinsAdj$HbA1c_Q <- qvalue(as.numeric(res_ProteinsAdj$HbA1c_pvalue))$qvalues
res_ProteinsAdj$WP_Q <- qvalue(as.numeric(res_ProteinsAdj$WP_pvalue))$qvalues

res_TargAdj <- res_Targ
res_TargAdj$Age_FDR <- p.adjust(as.numeric(res_TargAdj$Age_pvalue), method = "BH")
res_TargAdj$Sex_FDR <- p.adjust(as.numeric(res_TargAdj$Sex_pvalue), method = "BH")
res_TargAdj$BMI_FDR <- p.adjust(as.numeric(res_TargAdj$BMI_pvalue), method = "BH")
res_TargAdj$HbA1c_FDR <- p.adjust(as.numeric(res_TargAdj$HbA1c_pvalue), method = "BH")
res_TargAdj$WP_FDR <- p.adjust(as.numeric(res_TargAdj$WP_pvalue), method = "BH")
res_TargAdj$Age_Q <- qvalue(as.numeric(res_TargAdj$Age_pvalue), lambda = seq(0, max(as.numeric(res_TargAdj$Age_pvalue)), 0.05))$qvalues
res_TargAdj$Sex_Q <- qvalue(as.numeric(res_TargAdj$Sex_pvalue), lambda = seq(0, max(as.numeric(res_TargAdj$Sex_pvalue)), 0.05))$qvalues
res_TargAdj$BMI_Q <- qvalue(as.numeric(res_TargAdj$BMI_pvalue))$qvalues
res_TargAdj$HbA1c_Q <- qvalue(as.numeric(res_TargAdj$HbA1c_pvalue), lambda = seq(0, max(as.numeric(res_TargAdj$HbA1c_pvalue)), 0.05))$qvalues
res_TargAdj$WP_Q <- qvalue(as.numeric(res_TargAdj$WP_pvalue))$qvalues
res_TargAdj$CommonName <- sapply(strsplit(res_TargAdj$CommonName, "\\n"), `[`,1 )

res_UntargAdj <- res_Untarg
res_UntargAdj$Age_FDR <- p.adjust(as.numeric(res_UntargAdj$Age_pvalue), method = "BH")
res_UntargAdj$Sex_FDR <- p.adjust(as.numeric(res_UntargAdj$Sex_pvalue), method = "BH")
res_UntargAdj$BMI_FDR <- p.adjust(as.numeric(res_UntargAdj$BMI_pvalue), method = "BH")
res_UntargAdj$HbA1c_FDR <- p.adjust(as.numeric(res_UntargAdj$HbA1c_pvalue), method = "BH")
res_UntargAdj$WP_FDR <- p.adjust(as.numeric(res_UntargAdj$WP_pvalue), method = "BH")
res_UntargAdj$Age_Q <- qvalue(as.numeric(res_UntargAdj$Age_pvalue), lambda = seq(0, max(as.numeric(res_UntargAdj$Age_pvalue)), 0.05))$qvalues
res_UntargAdj$Sex_Q <- qvalue(as.numeric(res_UntargAdj$Sex_pvalue))$qvalues
res_UntargAdj$BMI_Q <- qvalue(as.numeric(res_UntargAdj$BMI_pvalue))$qvalues
res_UntargAdj$HbA1c_Q <- qvalue(as.numeric(res_UntargAdj$HbA1c_pvalue))$qvalues
res_UntargAdj$WP_Q <- qvalue(as.numeric(res_UntargAdj$WP_pvalue))$qvalues