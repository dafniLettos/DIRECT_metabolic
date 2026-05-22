# THIS FILE CONTAINS THE R CODE USED TO PLOT
# HISTOGRAMS OF P-VALUE DISTRIBUTIONS FOR AGE, SEX, BMI AND HbA1c ASSOCIATIONS WITH MOLECULAR PHENOTYPES
# VOLCANO PLOTS FOR THE ASSOCIATIONS
# ASSOCIATIONS' OVERLAP ACROSS TRAITS (AGE, SEX, BMI, HbA1c) - UpSETR
# BARPLOT FOR % ASSOCIATIONS
#-------------------------------------------------------------------	LIBRARIES/PACKAGES USED	-------------------------------------------------------------------------
library(tibble)
library(dplyr)
library(ggplot2)
library(reshape2)
library(UpSetR)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# DEFINE FUNCTIONS
`%!in%` <- Negate(`%in%`)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#		HISTOGRAMS OF NOMINAL P-VALUE DISTRIBUTIONS
data_plot <- data.frame(
    pVal = c(as.numeric(as.character(res_Genes$Age_pvalue)),as.numeric(as.character(res_Genes$BMI_pvalue)),as.numeric(as.character(res_Genes$Sex_pvalue)),as.numeric(as.character(res_Genes$HbA1c_pvalue)),as.numeric(as.character(res_Genes$WP_pvalue))),
    cov = c(rep("Age",nrow(res_Genes)),rep("BMI",nrow(res_Genes)),rep("Sex (M)",nrow(res_Genes)),rep("Fasting HbA1c",nrow(res_Genes)), rep("T2D Status",nrow(res_Genes))))
gplot(data = data_plot, aes(pVal)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "turquoise4") +
    scale_x_continuous(name = "p-value", labels = seq(0, 1, by = 0.25)) +
	labs(title="Nominal P-value Distributions for Gene Associations", x ="P-value", y = "") +
    facet_wrap(~cov) +
	theme( title = element_text( size = 14),
		axis.text = element_text( size = 10),
		axis.title = element_text( size = 12),
		strip.text = element_text(size = 12, face = "bold" ))
		
data_plot <- data.frame(
    pVal = c(as.numeric(as.character(res_Proteins$Age_pvalue)),as.numeric(as.character(res_Proteins$BMI_pvalue)),as.numeric(as.character(res_Proteins$Sex_pvalue)),as.numeric(as.character(res_Proteins$HbA1c_pvalue)),as.numeric(as.character(res_Proteins$WP_pvalue))),
    cov = c(rep("Age",nrow(res_Proteins)),rep("BMI",nrow(res_Proteins)),rep("Sex (M)",nrow(res_Proteins)),rep("Fasting HbA1c",nrow(res_Proteins)),rep("T2D Status",nrow(res_Proteins))))
ggplot(data = data_plot, aes(pVal)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "turquoise4") +
    scale_x_continuous(name = "p-value", labels = seq(0, 1, by = 0.25)) +
	labs(title="Nominal P-value Distributions for Protein Associations", x ="P-value", y = "") +
    facet_wrap(~cov) +
	theme( title = element_text( size = 14),
		axis.text = element_text( size = 10),
		axis.title = element_text( size = 12),
		strip.text = element_text(size = 12, face = "bold" ))

data_plot <- data.frame(
    pVal = c(as.numeric(as.character(res_Targ$Age_pvalue)),as.numeric(as.character(res_Targ$BMI_pvalue)),as.numeric(as.character(res_Targ$Sex_pvalue)),as.numeric(as.character(res_Targ$HbA1c_pvalue)),as.numeric(as.character(res_Targ$WP_pvalue))),
    cov = c(rep("Age",nrow(res_Targ)),rep("BMI",nrow(res_Targ)),rep("Sex (M)",nrow(res_Targ)),rep("Fasting HbA1c",nrow(res_Targ)),rep("T2D Status",nrow(res_Targ))))
ggplot(data = data_plot, aes(pVal)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "turquoise4") +
    scale_x_continuous(name = "p-value", labels = seq(0, 1, by = 0.25)) +
	labs(title="Nominal P-value Distributions for Targeted Metabolites Associations", x ="P-value", y = "") +
    facet_wrap(~cov) +
	theme( title = element_text( size = 14),
		axis.text = element_text( size = 10),
		axis.title = element_text( size = 12),
		strip.text = element_text(size = 12, face = "bold" ))

data_plot <- data.frame(
    pVal = c(as.numeric(as.character(res_Untarg$Age_pvalue)),as.numeric(as.character(res_Untarg$BMI_pvalue)),as.numeric(as.character(res_Untarg$Sex_pvalue)),as.numeric(as.character(res_Untarg$HbA1c_pvalue)),as.numeric(as.character(res_Untarg$WP_pvalue))),
    cov = c(rep("Age",nrow(res_Untarg)),rep("BMI",nrow(res_Untarg)),rep("Sex (M)",nrow(res_Untarg)),rep("Fasting HbA1c",nrow(res_Untarg)),rep("T2D Status",nrow(res_Untarg))))
ggplot(data = data_plot, aes(pVal)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), fill = "turquoise4") +
    scale_x_continuous(name = "p-value", labels = seq(0, 1, by = 0.25)) +
	labs(title="Nominal P-value Distributions for Untargeted Metabolites Associations", x ="P-value", y = "") +
    facet_wrap(~cov) +
	theme( title = element_text( size = 14),
		axis.text = element_text( size = 10),
		axis.title = element_text( size = 12),
		strip.text = element_text(size = 12, face = "bold" ))


#		VOLCANO PLOTS
res_GenesAdj$AgeSignifQ <- "NO"
res_GenesAdj$AgeSignifQ[as.double(res_GenesAdj$Age_Q) < 0.05] <- "SIGNIF"
res_GenesAdj$SexSignifQ <- "NO"
res_GenesAdj$SexSignifQ[as.double(res_GenesAdj$Sex_Q) < 0.05] <- "SIGNIF"
res_GenesAdj$BMISignifQ <- "NO"
res_GenesAdj$BMISignifQ[as.double(res_GenesAdj$BMI_Q) < 0.05] <- "SIGNIF"
res_GenesAdj$HbA1cSignifQ <- "NO"
res_GenesAdj$HbA1cSignifQ[as.double(res_GenesAdj$HbA1c_Q) < 0.05] <- "SIGNIF"
res_GenesAdj$AgeBetas <- "NO"
res_GenesAdj$AgeBetas[(as.double(res_GenesAdj$Age_Q) < 0.05)&(as.double(res_GenesAdj$Age) < 0)] <- "NEG"
res_GenesAdj$AgeBetas[(as.double(res_GenesAdj$Age_Q) < 0.05)&(as.double(res_GenesAdj$Age) > 0)] <- "POS"
res_GenesAdj$SexBetas <- "NO"
res_GenesAdj$SexBetas[(as.double(res_GenesAdj$Sex_Q) < 0.05)&(as.double(res_GenesAdj$Sex) < 0)] <- "NEG"
res_GenesAdj$SexBetas[(as.double(res_GenesAdj$Sex_Q) < 0.05)&(as.double(res_GenesAdj$Sex) > 0)] <- "POS"
res_GenesAdj$BMIBetas <- "NO"
res_GenesAdj$BMIBetas[(as.double(res_GenesAdj$BMI_Q) < 0.05)&(as.double(res_GenesAdj$BMI) < 0)] <- "NEG"
res_GenesAdj$BMIBetas[(as.double(res_GenesAdj$BMI_Q) < 0.05)&(as.double(res_GenesAdj$BMI) > 0)] <- "POS"
res_GenesAdj$HbA1cBetas <- "NO"
res_GenesAdj$HbA1cBetas[(as.double(res_GenesAdj$HbA1c_Q) < 0.05)&(as.double(res_GenesAdj$HbA1c) < 0)] <- "NEG"
res_GenesAdj$HbA1cBetas[(as.double(res_GenesAdj$HbA1c_Q) < 0.05)&(as.double(res_GenesAdj$HbA1c) > 0)] <- "POS"

ggplot(data=res_GenesAdj,aes(x=as.numeric(Age),y=-log10(as.numeric(Age_Q)),col=AgeBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="Age Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_GenesAdj,aes(x=as.numeric(Sex),y=-log10(as.numeric(Sex_Q)),col=SexBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="Sex (M) Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_GenesAdj,aes(x=as.numeric(BMI),y=-log10(as.numeric(BMI_Q)),col=BMIBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="BMI Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_GenesAdj,aes(x=as.numeric(HbA1c),y=-log10(as.numeric(HbA1c_Q)),col=HbA1cBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="HbA1c Estimates", x ="Betas", y = "-log10(P-value)")
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
res_ProteinsAdj$AgeSignifQ <- "NO"
res_ProteinsAdj$AgeSignifQ[as.double(res_ProteinsAdj$Age_Q) < 0.05] <- "SIGNIF"
res_ProteinsAdj$SexSignifQ <- "NO"
res_ProteinsAdj$SexSignifQ[as.double(res_ProteinsAdj$Sex_Q) < 0.05] <- "SIGNIF"
res_ProteinsAdj$BMISignifQ <- "NO"
res_ProteinsAdj$BMISignifQ[as.double(res_ProteinsAdj$BMI_Q) < 0.05] <- "SIGNIF"
res_ProteinsAdj$HbA1cSignifQ <- "NO"
res_ProteinsAdj$HbA1cSignifQ[as.double(res_ProteinsAdj$HbA1c_Q) < 0.05] <- "SIGNIF"
res_ProteinsAdj$AgeBetas <- "NO"
res_ProteinsAdj$AgeBetas[(as.double(res_ProteinsAdj$Age_Q) < 0.05)&(as.double(res_ProteinsAdj$Age) < 0)] <- "NEG"
res_ProteinsAdj$AgeBetas[(as.double(res_ProteinsAdj$Age_Q) < 0.05)&(as.double(res_ProteinsAdj$Age) > 0)] <- "POS"
res_ProteinsAdj$SexBetas <- "NO"
res_ProteinsAdj$SexBetas[(as.double(res_ProteinsAdj$Sex_Q) < 0.05)&(as.double(res_ProteinsAdj$Sex) < 0)] <- "NEG"
res_ProteinsAdj$SexBetas[(as.double(res_ProteinsAdj$Sex_Q) < 0.05)&(as.double(res_ProteinsAdj$Sex) > 0)] <- "POS"
res_ProteinsAdj$BMIBetas <- "NO"
res_ProteinsAdj$BMIBetas[(as.double(res_ProteinsAdj$BMI_Q) < 0.05)&(as.double(res_ProteinsAdj$BMI) < 0)] <- "NEG"
res_ProteinsAdj$BMIBetas[(as.double(res_ProteinsAdj$BMI_Q) < 0.05)&(as.double(res_ProteinsAdj$BMI) > 0)] <- "POS"
res_ProteinsAdj$HbA1cBetas <- "NO"
res_ProteinsAdj$HbA1cBetas[(as.double(res_ProteinsAdj$HbA1c_Q) < 0.05)&(as.double(res_ProteinsAdj$HbA1c) < 0)] <- "NEG"
res_ProteinsAdj$HbA1cBetas[(as.double(res_ProteinsAdj$HbA1c_Q) < 0.05)&(as.double(res_ProteinsAdj$HbA1c) > 0)] <- "POS"

ggplot(data=res_ProteinsAdj,aes(x=as.numeric(Age),y=-log10(as.numeric(Age_Q)),col=AgeBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="Age Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_ProteinsAdj,aes(x=as.numeric(Sex),y=-log10(as.numeric(Sex_Q)),col=SexBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="Sex (M) Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_ProteinsAdj,aes(x=as.numeric(BMI),y=-log10(as.numeric(BMI_Q)),col=BMIBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="BMI Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_ProteinsAdj,aes(x=as.numeric(HbA1c),y=-log10(as.numeric(HbA1c_Q)),col=HbA1cBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="HbA1c Estimates", x ="Betas", y = "-log10(P-value)")
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
res_TargAdj$AgeSignifQ <- "NO"
res_TargAdj$AgeSignifQ[as.double(res_TargAdj$Age_Q) < 0.05] <- "SIGNIF"
res_TargAdj$SexSignifQ <- "NO"
res_TargAdj$SexSignifQ[as.double(res_TargAdj$Sex_Q) < 0.05] <- "SIGNIF"
res_TargAdj$BMISignifQ <- "NO"
res_TargAdj$BMISignifQ[as.double(res_TargAdj$BMI_Q) < 0.05] <- "SIGNIF"
res_TargAdj$HbA1cSignifQ <- "NO"
res_TargAdj$HbA1cSignifQ[as.double(res_TargAdj$HbA1c_Q) < 0.05] <- "SIGNIF"
res_TargAdj$AgeBetas <- "NO"
res_TargAdj$AgeBetas[(as.double(res_TargAdj$Age_Q) < 0.05)&(as.double(res_TargAdj$Age) < 0)] <- "NEG"
res_TargAdj$AgeBetas[(as.double(res_TargAdj$Age_Q) < 0.05)&(as.double(res_TargAdj$Age) > 0)] <- "POS"
res_TargAdj$SexBetas <- "NO"
res_TargAdj$SexBetas[(as.double(res_TargAdj$Sex_Q) < 0.05)&(as.double(res_TargAdj$Sex) < 0)] <- "NEG"
res_TargAdj$SexBetas[(as.double(res_TargAdj$Sex_Q) < 0.05)&(as.double(res_TargAdj$Sex) > 0)] <- "POS"
res_TargAdj$BMIBetas <- "NO"
res_TargAdj$BMIBetas[(as.double(res_TargAdj$BMI_Q) < 0.05)&(as.double(res_TargAdj$BMI) < 0)] <- "NEG"
res_TargAdj$BMIBetas[(as.double(res_TargAdj$BMI_Q) < 0.05)&(as.double(res_TargAdj$BMI) > 0)] <- "POS"
res_TargAdj$HbA1cBetas <- "NO"
res_TargAdj$HbA1cBetas[(as.double(res_TargAdj$HbA1c_Q) < 0.05)&(as.double(res_TargAdj$HbA1c) < 0)] <- "NEG"
res_TargAdj$HbA1cBetas[(as.double(res_TargAdj$HbA1c_Q) < 0.05)&(as.double(res_TargAdj$HbA1c) > 0)] <- "POS"

ggplot(data=res_TargAdj,aes(x=as.numeric(Age),y=-log10(as.numeric(Age_Q)),col=AgeBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="Age Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_TargAdj,aes(x=as.numeric(Sex),y=-log10(as.numeric(Sex_Q)),col=SexBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="Sex (M) Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_TargAdj,aes(x=as.numeric(BMI),y=-log10(as.numeric(BMI_Q)),col=BMIBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="BMI Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_TargAdj,aes(x=as.numeric(HbA1c),y=-log10(as.numeric(HbA1c_Q)),col=HbA1cBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="HbA1c Estimates", x ="Betas", y = "-log10(P-value)")
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
res_UntargAdj$AgeSignifQ <- "NO"
res_UntargAdj$AgeSignifQ[as.double(res_UntargAdj$Age_Q) < 0.05] <- "SIGNIF"
res_UntargAdj$SexSignifQ <- "NO"
res_UntargAdj$SexSignifQ[as.double(res_UntargAdj$Sex_Q) < 0.05] <- "SIGNIF"
res_UntargAdj$BMISignifQ <- "NO"
res_UntargAdj$BMISignifQ[as.double(res_UntargAdj$BMI_Q) < 0.05] <- "SIGNIF"
res_UntargAdj$HbA1cSignifQ <- "NO"
res_UntargAdj$HbA1cSignifQ[as.double(res_UntargAdj$HbA1c_Q) < 0.05] <- "SIGNIF"
res_UntargAdj$AgeBetas <- "NO"
res_UntargAdj$AgeBetas[(as.double(res_UntargAdj$Age_Q) < 0.05)&(as.double(res_UntargAdj$Age) < 0)] <- "NEG"
res_UntargAdj$AgeBetas[(as.double(res_UntargAdj$Age_Q) < 0.05)&(as.double(res_UntargAdj$Age) > 0)] <- "POS"
res_UntargAdj$SexBetas <- "NO"
res_UntargAdj$SexBetas[(as.double(res_UntargAdj$Sex_Q) < 0.05)&(as.double(res_UntargAdj$Sex) < 0)] <- "NEG"
res_UntargAdj$SexBetas[(as.double(res_UntargAdj$Sex_Q) < 0.05)&(as.double(res_UntargAdj$Sex) > 0)] <- "POS"
res_UntargAdj$BMIBetas <- "NO"
res_UntargAdj$BMIBetas[(as.double(res_UntargAdj$BMI_Q) < 0.05)&(as.double(res_UntargAdj$BMI) < 0)] <- "NEG"
res_UntargAdj$BMIBetas[(as.double(res_UntargAdj$BMI_Q) < 0.05)&(as.double(res_UntargAdj$BMI) > 0)] <- "POS"
res_UntargAdj$HbA1cBetas <- "NO"
res_UntargAdj$HbA1cBetas[(as.double(res_UntargAdj$HbA1c_Q) < 0.05)&(as.double(res_UntargAdj$HbA1c) < 0)] <- "NEG"
res_UntargAdj$HbA1cBetas[(as.double(res_UntargAdj$HbA1c_Q) < 0.05)&(as.double(res_UntargAdj$HbA1c) > 0)] <- "POS"

ggplot(data=res_UntargAdj,aes(x=as.numeric(Age),y=-log10(as.numeric(Age_Q)),col=AgeBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="Age Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_UntargAdj,aes(x=as.numeric(Sex),y=-log10(as.numeric(Sex_Q)),col=SexBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="Sex (M) Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_UntargAdj,aes(x=as.numeric(BMI),y=-log10(as.numeric(BMI_Q)),col=BMIBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="BMI Estimates", x ="Betas", y = "-log10(P-value)")
ggplot(data=res_UntargAdj,aes(x=as.numeric(HbA1c),y=-log10(as.numeric(HbA1c_Q)),col=HbA1cBetas)) + geom_point() + theme_minimal() +
	geom_vline(xintercept=c(0,0), col="black") +
	scale_color_manual(values=setNames(c("maroon4", "gray70", "turquoise4"),c("POS","NO","NEG"))) +
	labs(title="HbA1c Estimates", x ="Betas", y = "-log10(P-value)")



#		UpSETR [Associations Overlap]
##	Gene Expression
lmGenesPlot <- data.frame(res_GenesHbA1cAdj$GeneID,res_GenesHbA1cAdj$Age_Q, res_GenesHbA1cAdj$BMI_Q, res_GenesHbA1cAdj$Sex_Q, res_GenesHbA1cAdj$HbA1c_Q)
colnames(lmGenesPlot) <- c("GeneID","Age_Q","BMI_Q","Sex_Q","HbA1c_Q")
lmGenesPlot$Age <- 0
lmGenesPlot$Age[lmGenesPlot$Age_Q < 0.05] <- 1
lmGenesPlot$BMI <- 0
lmGenesPlot$BMI[lmGenesPlot$BMI_Q < 0.05] <- 1
lmGenesPlot$Sex <- 0
lmGenesPlot$Sex[lmGenesPlot$Sex_Q < 0.05] <- 1
lmGenesPlot$HbA1c <- 0
lmGenesPlot$HbA1c[lmGenesPlot$HbA1c_Q < 0.05] <- 1
lmGenesPlot[,c(2:5)] <- NULL

upset(lmGenesPlot, sets = c("Age", "BMI", "Sex", "HbA1c"), main.bar.color="turquoise4" , matrix.color="turquoise4", sets.bar.color = "maroon4",
keep.order = TRUE, mainbar.y.label = "Common Genes", sets.x.label = "Significant Gene Count", 
text.scale = c(1.7,1.3,1.3,1.3,1.4,1.5), order.by = "degree", empty.intersections = "on")
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##	Proteins
lmProteinsPlot <- data.frame(res_ProteinsHbA1cAdj$GeneID,res_ProteinsHbA1cAdj$Age_Q, res_ProteinsHbA1cAdj$BMI_Q, res_ProteinsHbA1cAdj$Sex_Q, res_ProteinsHbA1cAdj$HbA1c_Q)
colnames(lmProteinsPlot) <- c("ProteinID","Age_Q","BMI_Q","Sex_Q","HbA1c_Q")
lmProteinsPlot$Age <- 0
lmProteinsPlot$Age[lmProteinsPlot$Age_Q < 0.05] <- 1
lmProteinsPlot$BMI <- 0
lmProteinsPlot$BMI[lmProteinsPlot$BMI_Q < 0.05] <- 1
lmProteinsPlot$Sex <- 0
lmProteinsPlot$Sex[lmProteinsPlot$Sex_Q < 0.05] <- 1
lmProteinsPlot$HbA1c <- 0
lmProteinsPlot$HbA1c[lmProteinsPlot$HbA1c_Q < 0.05] <- 1
lmProteinsPlot[,c(2:5)] <- NULL

upset(lmProteinsPlot, sets = c("Age", "BMI", "Sex", "HbA1c"), main.bar.color="turquoise4" , matrix.color="turquoise4", sets.bar.color = "maroon4",
keep.order = TRUE, mainbar.y.label = "Common Proteins", sets.x.label = "Significant Protein Count", 
text.scale = c(1.7,1.3,1.3,1.3,1.4,1.5), order.by = "degree", empty.intersections = "on")
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
##	Metabolites
res_UntargAdj_tmp <- res_UntargAdj[,c(1,3,6:40)]
colnames(res_UntargAdj_tmp)[1:2] <- c("Metabolite","CommonName")
res_Metabolites_merged <- rbind(res_TargAdj[,c(1:2,4:38)],res_UntargAdj_tmp)
lmMetabPlot <- data.frame(res_Metabolites_merged$Metabolite,res_Metabolites_merged$Age_Q, res_Metabolites_merged$BMI_Q, res_Metabolites_merged$Sex_Q, res_Metabolites_merged$HbA1c_Q)
colnames(lmMetabPlot) <- c("Metabolite","Age_Q","BMI_Q","Sex_Q","HbA1c_Q")
lmMetabPlot$Age <- 0
lmMetabPlot$Age[lmMetabPlot$Age_Q < 0.05] <- 1
lmMetabPlot$BMI <- 0
lmMetabPlot$BMI[lmMetabPlot$BMI_Q < 0.05] <- 1
lmMetabPlot$Sex <- 0
lmMetabPlot$Sex[lmMetabPlot$Sex_Q < 0.05] <- 1
lmMetabPlot$HbA1c <- 0
lmMetabPlot$HbA1c[lmMetabPlot$HbA1c_Q < 0.05] <- 1
lmMetabPlot[,c(2:5)] <- NULL

upset(lmMetabPlot, sets = c("Age", "BMI", "Sex", "HbA1c"), main.bar.color="turquoise4" , matrix.color="turquoise4", sets.bar.color = "maroon4",
keep.order = TRUE, mainbar.y.label = "Targeted Metabolites", sets.x.label = "Significant Metabolite Count", 
text.scale = c(1.7,1.3,1.3,1.3,1.4,1.5), order.by = "degree", empty.intersections = "on")
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#		BARPLOT
phenotype_counts <- matrix(0,ncol=4,nrow=5)
colnames(phenotype_counts) <- c("Age","Sex","BMI","HbA1c")
rownames(phenotype_counts) <- c("Gene Expression","Proteins","Metab_Targ","Metab_Untarg","Metab_Merged")
phenotype_counts[1,1] <- nrow(subset(res_GenesAdj, Age_Q < 0.05))
phenotype_counts[1,2] <- nrow(subset(res_GenesAdj, Sex_Q < 0.05))
phenotype_counts[1,3] <- nrow(subset(res_GenesAdj, BMI_Q < 0.05))
phenotype_counts[1,4] <- nrow(subset(res_GenesAdj, HbA1c_Q < 0.05))
phenotype_counts[2,1] <- nrow(subset(res_ProteinsAdj, Age_Q < 0.05))
phenotype_counts[2,2] <- nrow(subset(res_ProteinsAdj, Sex_Q < 0.05))
phenotype_counts[2,3] <- nrow(subset(res_ProteinsAdj, BMI_Q < 0.05))
phenotype_counts[2,4] <- nrow(subset(res_ProteinsAdj, HbA1c_Q < 0.05))
phenotype_counts[3,1] <- nrow(subset(res_TargAdj, Age_Q < 0.05))
phenotype_counts[3,2] <- nrow(subset(res_TargAdj, Sex_Q < 0.05))
phenotype_counts[3,3] <- nrow(subset(res_TargAdj, BMI_Q < 0.05))
phenotype_counts[3,4] <- nrow(subset(res_TargAdj, HbA1c_Q < 0.05))
phenotype_counts[4,1] <- nrow(subset(res_UntargAdj, Age_Q < 0.05))
phenotype_counts[4,2] <- nrow(subset(res_UntargAdj, Sex_Q < 0.05))
phenotype_counts[4,3] <- nrow(subset(res_UntargAdj, BMI_Q < 0.05))
phenotype_counts[4,4] <- nrow(subset(res_UntargAdj, HbA1c_Q < 0.05))
phenotype_counts[5,1] <- length(subset(res_TargAdj, select=(Metabolite), subset=(Age_Q < 0.05))[,1]) +  length(subset(res_UntargAdj, select=(COMP_ID), subset=(Age_Q < 0.05))[,1])
phenotype_counts[5,2] <- length(subset(res_TargAdj, select=(Metabolite), subset=(Sex_Q < 0.05))[,1]) +  length(subset(res_UntargAdj, select=(COMP_ID), subset=(Sex_Q < 0.05))[,1])
phenotype_counts[5,3] <- length(subset(res_TargAdj, select=(Metabolite), subset=(BMI_Q < 0.05))[,1]) +  length(subset(res_UntargAdj, select=(COMP_ID), subset=(BMI_Q < 0.05))[,1])
phenotype_counts[5,4] <- length(subset(res_TargAdj, select=(Metabolite), subset=(HbA1c_Q < 0.05))[,1]) +  length(subset(res_UntargAdj, select=(COMP_ID), subset=(HbA1c_Q < 0.05))[,1])

phenotype_percentages <- matrix(0,ncol=4,nrow=5)
colnames(phenotype_percentages) <- c("Age","Sex","BMI","HbA1c")
rownames(phenotype_percentages) <- c("Gene Expression","Proteins","Metab_T","Metab_U","Metab_Merged")
phenotype_percentages[1,1] <- phenotype_counts[1,1]*100/16209
phenotype_percentages[1,2] <- phenotype_counts[1,2]*100/16209
phenotype_percentages[1,3] <- phenotype_counts[1,3]*100/16209
phenotype_percentages[1,4] <- phenotype_counts[1,4]*100/16209
phenotype_percentages[2,1] <- phenotype_counts[2,1]*100/373
phenotype_percentages[2,2] <- phenotype_counts[2,2]*100/373
phenotype_percentages[2,3] <- phenotype_counts[2,3]*100/373
phenotype_percentages[2,4] <- phenotype_counts[2,4]*100/373
phenotype_percentages[3,1] <- phenotype_counts[3,1]*100/116
phenotype_percentages[3,2] <- phenotype_counts[3,2]*100/116
phenotype_percentages[3,3] <- phenotype_counts[3,3]*100/116
phenotype_percentages[3,4] <- phenotype_counts[3,4]*100/116
phenotype_percentages[4,1] <- phenotype_counts[4,1]*100/233
phenotype_percentages[4,2] <- phenotype_counts[4,2]*100/233
phenotype_percentages[4,3] <- phenotype_counts[4,3]*100/233
phenotype_percentages[4,4] <- phenotype_counts[4,4]*100/233
phenotype_percentages[5,1] <- phenotype_counts[5,1]*100/349
phenotype_percentages[5,2] <- phenotype_counts[5,2]*100/349
phenotype_percentages[5,3] <- phenotype_counts[5,3]*100/349
phenotype_percentages[5,4] <- phenotype_counts[5,4]*100/349
phenotype_percentages2 <- data.frame(Associations = c(phenotype_percentages[1:5,1],phenotype_percentages[1:5,2],phenotype_percentages[1:5,3],phenotype_percentages[1:5,4]), 
												covariate = rep(c("Age", "Sex", "BMI","HbA1c"), each = 5),
												phenotype = c("Gene Expression","Proteins","Metabolites (Targeted)","Metabolites (Untargeted)","Metabolites (Merged)"))
  
ggplot(phenotype_percentages2, aes(x = factor(covariate, levels = c("Age", "Sex", "BMI","HbA1c")), y =Associations, fill = factor(phenotype, levels = c("Gene Expression","Proteins","Metabolites (Targeted)","Metabolites (Untargeted)","Metabolites (Merged)")))) +
	geom_bar(stat = "identity", position = "dodge") +
	labs(x=NULL, y = "% Significant Associations", fill = NULL ) +
    scale_fill_manual(values=c("#CC6666", "#9999CC", "palegreen", "palegreen3", "palegreen4")) +
	theme( title = element_text( size = 16),
		axis.text = element_text( size = 14, face = "bold" ),
		axis.title = element_text( size = 16),
		legend.text = element_text(size = 12))