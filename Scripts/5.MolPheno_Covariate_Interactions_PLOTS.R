# THIS FILE CONTAINS THE R CODE USED TO PLOT
# HISTOGRAMS OF P-VALUE DISTRIBUTIONS FOR AGE, SEX, BMI AND HbA1c ASSOCIATIONS WITH MOLECULAR PHENOTYPES
# VOLCANO PLOTS FOR THE ASSOCIATIONS
# ASSOCIATIONS' OVERLAP ACROSS TRAITS (AGE, SEX, BMI, HbA1c) - UpSETR
# BARPLOT FOR % ASSOCIATIONS
#-------------------------------------------------------------------	LIBRARIES/PACKAGES USED	-------------------------------------------------------------------------
library(tibble)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(UpSetR)
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lm_genes <- read.delim("lmResultsAdjGenes_HbA1cm.txt", header = TRUE)
lm_proteins <- read.delim("lmResultsAdjProteins_HbA1cm.txt", header = TRUE)
lm_metabT <- read.delim("lmResultsAdjTarMet_HbA1cm.txt", header = TRUE)
lm_metabU <- read.delim("lmResultsAdjUntarMet_HbA1cm.txt", header = TRUE)

lm_genes$GeneID <- trimws(lm_genes$GeneID)
lm_genes$GeneName <- trimws(lm_genes$GeneName)
lm_genes_age <- subset(lm_genes, subset=(Age_Q < 0.05))
lm_proteins$ProteinID <- trimws(lm_proteins$ProteinID)
lm_proteins$GeneID2 <- trimws(lm_proteins$GeneID2)
lm_proteins$GeneName <- trimws(lm_proteins$GeneName)
lm_proteins_age <- subset(lm_proteins, subset=(Age_Q < 0.05))
lm_metabT$CommonName <- trimws(lm_metabT$CommonName)
lm_metabU$BIOCHEMICAL <- trimws(lm_metabU$BIOCHEMICAL)
lm_metabU$BIOCHEMICAL2 <- trimws(lm_metabU$BIOCHEMICAL2)
lm_metabT_age <- subset(lm_metabT, subset=(Age_Q < 0.05))
lm_metabU_age <- subset(lm_metabU, subset=(Age_Q < 0.05))
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#		SCATTERPLOTS - AGE*SEX
##	Gene Expression
ggplot(lm_genes_age,aes(x=as.numeric(Age),y=as.numeric(Sex),col=SexSignifQ, label=GeneName)) + geom_point() +
	geom_vline(xintercept=c(0,0), col="black") +
	geom_hline(yintercept=c(0,0), col="black") +
	scale_color_manual(values=c("gray70", "darkred")) +
	guides(col = FALSE) +
	geom_label_repel(data = subset(lm_genes_age, subset=(Age_Q < 0.01 & Sex_Q < 0.01)),
		size = 2, max.overlaps = 61,
        box.padding   = 0.01, 
        point.padding = 0.01,
        segment.color = 'grey50') +
	labs(title="Age Associations in Male & Female Individuals - Gene Expression", x ="Age Beta", y = "Sex Beta") +
		theme(title = element_text( size = 14), panel.grid = element_blank(),
		axis.text = element_text( size = 12),
		axis.title = element_text( size = 14))
dim(subset(lm_genes_age, subset=(Age_Q < 0.01 & Sex_Q < 0.01)))
#[1] 4580   41

sum((lm_genes_age$Age > 0)&(lm_genes_age$Sex > 0))	3194
sum((lm_genes_age$Age > 0)&(lm_genes_age$Sex < 0))	747
sum((lm_genes_age$Age < 0)&(lm_genes_age$Sex > 0))	2368
sum((lm_genes_age$Age < 0)&(lm_genes_age$Sex < 0))	4334
sum((lm_genes_age$Age > 0)&(lm_genes_age$Sex > 0)&(lm_genes_age$SexSignifQ == " SIGNIF"))	2725
sum((lm_genes_age$Age > 0)&(lm_genes_age$Sex < 0)&(lm_genes_age$SexSignifQ == " SIGNIF"))	397
sum((lm_genes_age$Age < 0)&(lm_genes_age$Sex > 0)&(lm_genes_age$SexSignifQ == " SIGNIF"))	1075
sum((lm_genes_age$Age < 0)&(lm_genes_age$Sex < 0)&(lm_genes_age$SexSignifQ == " SIGNIF"))	2826

for(i in 1:nrow(lm_genes)){
	if(lm_genes$GeneID[i] %in% lm_genes_AGExSEX$GeneID){
		lm_genes$AGExSEX_Q[i] <- with(lm_genes_AGExSEX, Age.Sex_Q[GeneID == lm_genes$GeneID[i]])}
	else{lm_genes$AGExSEX_Q[i] <- 1}
}
lm_genes$AGExSEX_SignifQ <- "NO"
lm_genes$AGExSEX_SignifQ[lm_genes$AGExSEX_Q < 0.05] <- "SIGNIF"
ggplot(subset(lm_genes, subset=(AGExSEX_SignifQ == "SIGNIF")), aes(x=as.numeric(Age), y=as.numeric(Sex), col=AgeSignifQ, label=GeneName)) + geom_point() +
	geom_vline(xintercept=c(0,0), col="black") +
	geom_hline(yintercept=c(0,0), col="black") +
	scale_color_manual(labels = c("Not Significant", "Significant"), values=c("black", "darkred")) +
	geom_label_repel(data = subset(lm_genes, subset=(AGExSEX_Q < 0.025)),
		size = 2, max.overlaps = 1500,
        box.padding   = 0.01, 
        point.padding = 0.01,
        segment.color = 'grey50') +
	labs(title="Sex-specific Effects of Age - Gene Expression", x ="Age Beta", y = "Sex Beta", color = "AGE") +
		theme(title = element_text( size = 14), panel.grid = element_blank(),
		axis.text = element_text( size = 12),
		axis.title = element_text( size = 14))	


##	Proteins
ggplot(lm_proteins_age,aes(x=as.numeric(Age),y=as.numeric(Sex),col=SexSignifQ, label=GeneName)) + geom_point() +
	geom_vline(xintercept=c(0,0), col="black") +
	geom_hline(yintercept=c(0,0), col="black") +
	scale_color_manual(values=c("gray70", "mediumpurple4")) +
	guides(col = FALSE) +
	geom_label_repel(data = subset(lm_proteins_age, subset=(Age_Q < 0.01 & Sex_Q < 0.01)),
		size = 2, max.overlaps = 30,
        box.padding   = 0.01, 
        point.padding = 0.01,
        segment.color = 'grey50') +
	labs(title="Age Associations in Male & Female Individuals - Proteins", x ="Age Beta", y = "Sex Beta") +
		theme(title = element_text( size = 14), panel.grid = element_blank(),
		axis.text = element_text( size = 12),
		axis.title = element_text( size = 14))
dim(subset(lm_proteins_age, subset=(Age_Q < 0.01 & Sex_Q < 0.01)))
#[1] 149   41
	
sum((lm_proteins_age$Age > 0)&(lm_proteins_age$Sex > 0))	105
sum((lm_proteins_age$Age > 0)&(lm_proteins_age$Sex < 0))	110
sum((lm_proteins_age$Age < 0)&(lm_proteins_age$Sex > 0))	20
sum((lm_proteins_age$Age < 0)&(lm_proteins_age$Sex < 0))	32
sum((lm_proteins_age$Age > 0)&(lm_proteins_age$Sex > 0)&(lm_proteins_age$SexSignifQ == " SIGNIF"))	75
sum((lm_proteins_age$Age > 0)&(lm_proteins_age$Sex < 0)&(lm_proteins_age$SexSignifQ == " SIGNIF"))	82
sum((lm_proteins_age$Age < 0)&(lm_proteins_age$Sex > 0)&(lm_proteins_age$SexSignifQ == " SIGNIF"))	16
sum((lm_proteins_age$Age < 0)&(lm_proteins_age$Sex < 0)&(lm_proteins_age$SexSignifQ == " SIGNIF"))	25

for(i in 1:nrow(lm_proteins)){
	if(lm_proteins$ProteinID[i] %in% lm_proteins_AGExSEX$ProteinID){
		lm_proteins$AGExSEX_Q[i] <- with(lm_proteins_AGExSEX, Age.Sex_Q[ProteinID == lm_proteins$ProteinID[i]])}
	else{lm_proteins$AGExSEX_Q[i] <- 1}
}
lm_proteins$AGExSEX_SignifQ <- "NO"
lm_proteins$AGExSEX_SignifQ[lm_proteins$AGExSEX_Q < 0.05] <- "SIGNIF"
ggplot(subset(lm_proteins, subset=(AGExSEX_SignifQ == "SIGNIF")),aes(x=as.numeric(Age),y=as.numeric(Sex),col=AgeSignifQ, label=GeneName)) + geom_point() +
	geom_vline(xintercept=c(0,0), col="black") +
	geom_hline(yintercept=c(0,0), col="black") +
	scale_color_manual(labels = c("Not Significant", "Significant"), values=c("black", "mediumpurple4")) +
	geom_label_repel(data = subset(lm_proteins, subset=(AGExSEX_Q < 0.05)),
		size = 2, max.overlaps = 100,
        box.padding   = 0.01, 
        point.padding = 0.01,
        segment.color = 'grey50') +
	labs(title="Sex-specific Effects of Age - Proteins", x ="Age Beta", y = "Sex Beta", color = "AGE") +
		theme(title = element_text( size = 14), panel.grid = element_blank(),
		axis.text = element_text( size = 12),
		axis.title = element_text( size = 14))


##	Metabolites
lm_metabT_age_temp <- age_metabT_HMDB[,c(2:3,5:6,25:26,31)] #102
lm_metabU_age_temp <- age_metabU_HMDB[,c(3:4,7:8,27:28,33)] #138
colnames(lm_metabT_age_temp) <- c("Metabolite","HMDB","Age","Sex","Age_Q", "Sex_Q", "SexSignifQ")
colnames(lm_metabU_age_temp) <- c("Metabolite","HMDB","Age","Sex","Age_Q", "Sex_Q", "SexSignifQ")
lm_metab_age_merged <- rbind(lm_metabT_age_temp, lm_metabU_age_temp) #240

ggplot(lm_metab_age_merged,aes(x=as.numeric(Age),y=as.numeric(Sex),col=SexSignifQ, label=Metabolite)) + geom_point() +
	geom_vline(xintercept=c(0,0), col="black") +
	geom_hline(yintercept=c(0,0), col="black") +
	scale_color_manual(values=c("gray70", "darkgreen")) +
	guides(col = FALSE) +
	geom_label_repel(data = subset(lm_metab_age_merged, subset=(Age_Q < 0.01 & Sex_Q < 0.01)),
		size = 2, max.overlaps = 10,
		box.padding   = 0.01, 
        point.padding = 0.01,
        segment.color = 'grey50') +
	labs(title="Age Associations in Male & Female Individuals - Metabolites", x ="Age Beta", y = "Sex Beta") +
		theme(title = element_text( size = 14), panel.grid = element_blank(),
		axis.text = element_text( size = 12),
		axis.title = element_text( size = 14))	
dim(subset(lm_metab_age_merged, subset=(Age_Q < 0.01 & Sex_Q < 0.01)))
#[1] 149   7
	
sum((lm_metab_age_merged$Age > 0)&(lm_metab_age_merged$Sex > 0))	65
sum((lm_metab_age_merged$Age > 0)&(lm_metab_age_merged$Sex < 0))	98
sum((lm_metab_age_merged$Age < 0)&(lm_metab_age_merged$Sex > 0))	48
sum((lm_metab_age_merged$Age < 0)&(lm_metab_age_merged$Sex < 0))	29
sum((lm_metab_age_merged$Age > 0)&(lm_metab_age_merged$Sex > 0)&(lm_metab_age_merged == " SIGNIF"))	55
sum((lm_metab_age_merged$Age > 0)&(lm_metab_age_merged$Sex < 0)&(lm_metab_age_merged == " SIGNIF"))	80
sum((lm_metab_age_merged$Age < 0)&(lm_metab_age_merged$Sex > 0)&(lm_metab_age_merged == " SIGNIF"))	47
sum((lm_metab_age_merged$Age < 0)&(lm_metab_age_merged$Sex < 0)&(lm_metab_age_merged == " SIGNIF"))	28

for(i in 1:nrow(lm_metabT)){
	if(lm_metabT$Metabolite[i] %in% lm_metabT_AGExSEX$Metabolite){
		lm_metabT$AGExSEX_Q[i] <- with(lm_metabT_AGExSEX, Age.Sex_Q[Metabolite == lm_metabT$Metabolite[i]])}
	else{lm_metabT$AGExSEX_Q[i] <- 1}
}
lm_metabT$AGExSEX_SignifQ <- "NO"
lm_metabT$AGExSEX_SignifQ[lm_metabT$AGExSEX_Q < 0.05] <- "SIGNIF"
for(i in 1:nrow(lm_metabU)){
	if(lm_metabU$COMP_ID[i] %in% lm_metabU_AGExSEX$COMP_ID){
		lm_metabU$AGExSEX_Q[i] <- with(lm_metabU_AGExSEX, Age.Sex_Q[COMP_ID == lm_metabU$COMP_ID[i]])}
	else{lm_metabU$AGExSEX_Q[i] <- 1}
}
lm_metabU$AGExSEX_SignifQ <- "NO"
lm_metabU$AGExSEX_SignifQ[lm_metabU$AGExSEX_Q < 0.05] <- "SIGNIF"

metabT_ageSex_temp <- lm_metabT[,c(2,5:6,25:26,30:31,40:41)]
metabU_ageSex_temp <- lm_metabU[,c(3,7:8,27:28,32:33,42:43)]
colnames(metabT_ageSex_temp) <- c("Metabolite","Age","Sex","Age_Q", "Sex_Q", "AgeSignifQ", "SexSignifQ", "AGExSEX_Q", "AGExSEX_SignifQ")
colnames(metabU_ageSex_temp) <- c("Metabolite","Age","Sex","Age_Q", "Sex_Q", "AgeSignifQ", "SexSignifQ", "AGExSEX_Q", "AGExSEX_SignifQ")
metab_ageSex_merged <- rbind(metabT_ageSex_temp, metabU_ageSex_temp)	

ggplot(subset(metab_ageSex_merged, subset=(AGExSEX_SignifQ == "SIGNIF")),aes(x=as.numeric(Age),y=as.numeric(Sex),col=AgeSignifQ, label=Metabolite)) + geom_point() +
	geom_vline(xintercept=c(0,0), col="black") +
	geom_hline(yintercept=c(0,0), col="black") +
	scale_color_manual(labels = c("Not Significant", "Significant"), values=c("black", "darkgreen")) +
	geom_label_repel(data = subset(metab_ageSex_merged, subset=(AGExSEX_Q < 0.05)), 
		size = 2, max.overlaps = 100,
        box.padding   = 0.01, 
        point.padding = 0.01,
        segment.color = 'grey50') +
	labs(title="Sex-specific Effects of Age - Metabolites", x ="Age Beta", y = "Sex Beta", color="AGE") +
		theme(title = element_text( size = 14), panel.grid = element_blank(),
		axis.text = element_text( size = 12),
		axis.title = element_text( size = 14))
#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lm_genes_age_BMI <- subset(lm_genes, subset=(Age_Q < 0.05 & BMI_Q < 0.05))
lm_proteins_age_BMI <- subset(lm_proteins,  subset=(Age_Q < 0.05 & BMI_Q < 0.05))
lm_metabT_age_BMI <- subset(lm_metabT,  subset=(Age_Q < 0.05 & BMI_Q < 0.05))
lm_metabU_age_BMI <- subset(lm_metabU,  subset=(Age_Q < 0.05 & BMI_Q < 0.05))
lm_metab_merged_age_BMI <- subset(lm_metab_merged,  subset=(Age_Q < 0.05 & BMI_Q < 0.05))

lm_genes_age_HbA1c <- subset(lm_genes, subset=(Age_Q < 0.05 & HbA1c_Q < 0.05))
lm_proteins_age_HbA1c <- subset(lm_proteins,  subset=(Age_Q < 0.05 & HbA1c_Q < 0.05))
lm_metabT_age_HbA1c <- subset(lm_metabT,  subset=(Age_Q < 0.05 & HbA1c_Q < 0.05))
lm_metabU_age_HbA1c <- subset(lm_metabU,  subset=(Age_Q < 0.05 & HbA1c_Q < 0.05))
lm_metab_merged_age_HbA1c <- subset(lm_metab_merged,  subset=(Age_Q < 0.05 & HbA1c_Q < 0.05))


#		SCATTERPLOTS - AGE, BMI, HbA1c
##	Gene Expression
ggplot(lm_genes_age_BMI,aes(x=as.numeric(Age),y=as.numeric(BMI),col=HbA1cSignifQ, label=GeneName)) + geom_point() +
	geom_vline(xintercept=c(0,0), col="black") +
	geom_hline(yintercept=c(0,0), col="black") +
	scale_color_manual(values=c("gray70", "darkred")) +
	guides(col = FALSE) +
	geom_label_repel(data = subset(lm_genes_age_BMI, subset=(Age_Q < 0.01 & BMI_Q < 0.01 & HbA1c_Q < 0.01)),
		size = 2, max.overlaps = 20,
        box.padding   = 0.01, 
        point.padding = 0.01,
        segment.color = 'grey50') +
	labs(title="Age and BMI Associations: Gene Expression", x ="Age Beta", y = "BMI Beta") +
		theme(title = element_text( size = 14), panel.grid = element_blank(),
		axis.text = element_text( size = 12),
		axis.title = element_text( size = 14))

sum((lm_genes_age_BMI$Age > 0)&(lm_genes_age_BMI$BMI > 0))	# 964
sum((lm_genes_age_BMI$Age > 0)&(lm_genes_age_BMI$BMI < 0))	# 1153
sum((lm_genes_age_BMI$Age < 0)&(lm_genes_age_BMI$BMI > 0))	# 457
sum((lm_genes_age_BMI$Age < 0)&(lm_genes_age_BMI$BMI < 0))	# 3124
sum((lm_genes_age_BMI$Age > 0)&(lm_genes_age_BMI$BMI > 0)&(lm_genes_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 332
sum((lm_genes_age_BMI$Age > 0)&(lm_genes_age_BMI$BMI < 0)&(lm_genes_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 286
sum((lm_genes_age_BMI$Age < 0)&(lm_genes_age_BMI$BMI > 0)&(lm_genes_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 134
sum((lm_genes_age_BMI$Age < 0)&(lm_genes_age_BMI$BMI < 0)&(lm_genes_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 1091

ggplot(lm_genes_age_HbA1c,aes(x=as.numeric(Age),y=as.numeric(HbA1c),col=BMISignifQ, label=GeneName)) + geom_point() +
  geom_vline(xintercept=c(0,0), col="black") +
  geom_hline(yintercept=c(0,0), col="black") +
  scale_color_manual(values=c("gray70", "darkred")) +
  guides(col = FALSE) +
  geom_label_repel(data = subset(lm_genes_age_HbA1c, subset=(Age_Q < 0.01 & BMI_Q < 0.01 & HbA1c_Q < 0.01)),
                   size = 2, max.overlaps = 20,
                   box.padding   = 0.01, 
                   point.padding = 0.01,
                   segment.color = 'grey50') +
  labs(title="Age and HbA1c Associations: Gene Expression", x ="Age Beta", y = "HbA1c Beta") +
  theme(title = element_text( size = 14), panel.grid = element_blank(),
        axis.text = element_text( size = 12),
        axis.title = element_text( size = 14))

sum((lm_genes_age_HbA1c$Age > 0)&(lm_genes_age_HbA1c$HbA1c > 0))	# 667
sum((lm_genes_age_HbA1c$Age > 0)&(lm_genes_age_HbA1c$HbA1c < 0))	# 425
sum((lm_genes_age_HbA1c$Age < 0)&(lm_genes_age_HbA1c$HbA1c > 0))	# 147
sum((lm_genes_age_HbA1c$Age < 0)&(lm_genes_age_HbA1c$HbA1c < 0))	# 2120
sum((lm_genes_age_HbA1c$Age > 0)&(lm_genes_age_HbA1c$HbA1c > 0)&(lm_genes_age_HbA1c$BMISignifQ == " SIGNIF"))	# 399
sum((lm_genes_age_HbA1c$Age > 0)&(lm_genes_age_HbA1c$HbA1c < 0)&(lm_genes_age_HbA1c$BMISignifQ == " SIGNIF"))	# 219
sum((lm_genes_age_HbA1c$Age < 0)&(lm_genes_age_HbA1c$HbA1c > 0)&(lm_genes_age_HbA1c$BMISignifQ == " SIGNIF"))	# 86
sum((lm_genes_age_HbA1c$Age < 0)&(lm_genes_age_HbA1c$HbA1c < 0)&(lm_genes_age_HbA1c$BMISignifQ == " SIGNIF"))	# 1139



##	Proteins
prot_BMI_by_HbA1c <- c("FETUB","SRCRB4D","FABP4","PON3","LEP","IGFBP1","IGFBP2","ANG")

ggplot(lm_proteins_age_BMI,aes(x=as.numeric(Age),y=as.numeric(BMI),col=HbA1cSignifQ, label=GeneName)) + geom_point() +
  geom_vline(xintercept=c(0,0), col="black") +
  geom_hline(yintercept=c(0,0), col="black") +
  scale_color_manual(values=c("gray70", "mediumpurple4")) +
  guides(col = FALSE) +
  geom_label_repel(data = subset(lm_proteins_age_BMI, subset=(abs(Age)>0.01 & abs(BMI)>0.04 & HbA1c_Q < 0.05 | GeneName %in% prot_BMI_by_HbA1c)),
                   size = 2, max.overlaps = 20,
                   box.padding   = 0.01, 
                   point.padding = 0.01,
                   segment.color = 'grey50') +
  labs(title="Age and BMI Associations: Protein Levels", x ="Age Beta", y = "BMI Beta") +
  theme(title = element_text( size = 14), panel.grid = element_blank(),
        axis.text = element_text( size = 12),
        axis.title = element_text( size = 14))

sum((lm_proteins_age_BMI$Age > 0)&(lm_proteins_age_BMI$BMI > 0))	# 163
sum((lm_proteins_age_BMI$Age > 0)&(lm_proteins_age_BMI$BMI < 0))	# 46
sum((lm_proteins_age_BMI$Age < 0)&(lm_proteins_age_BMI$BMI > 0))	# 36
sum((lm_proteins_age_BMI$Age < 0)&(lm_proteins_age_BMI$BMI < 0))	# 15
sum((lm_proteins_age_BMI$Age > 0)&(lm_proteins_age_BMI$BMI > 0)&(lm_proteins_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 90
sum((lm_proteins_age_BMI$Age > 0)&(lm_proteins_age_BMI$BMI < 0)&(lm_proteins_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 21
sum((lm_proteins_age_BMI$Age < 0)&(lm_proteins_age_BMI$BMI > 0)&(lm_proteins_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 23
sum((lm_proteins_age_BMI$Age < 0)&(lm_proteins_age_BMI$BMI < 0)&(lm_proteins_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 6

ggplot(lm_proteins_age_HbA1c,aes(x=as.numeric(Age),y=as.numeric(HbA1c),col=BMISignifQ, label=GeneName)) + geom_point() +
  geom_vline(xintercept=c(0,0), col="black") +
  geom_hline(yintercept=c(0,0), col="black") +
  scale_color_manual(values=c("gray70", "mediumpurple4")) +
  guides(col = FALSE) +
  geom_label_repel(data = subset(lm_proteins_age_HbA1c, subset=(abs(Age)>0.01 & abs(BMI)>0.01 & BMI_Q < 0.05 | GeneName %in% prot_BMI_by_HbA1c)),
                   size = 2, max.overlaps = 20,
                   box.padding   = 0.01, 
                   point.padding = 0.01,
                   segment.color = 'grey50') +
  labs(title="Age and HbA1c Associations: Protein Levels", x ="Age Beta", y = "HbA1c Beta") +
  theme(title = element_text( size = 14), panel.grid = element_blank(),
        axis.text = element_text( size = 12),
        axis.title = element_text( size = 14))

sum((lm_proteins_age_HbA1c$Age > 0)&(lm_proteins_age_HbA1c$HbA1c > 0))	# 99
sum((lm_proteins_age_HbA1c$Age > 0)&(lm_proteins_age_HbA1c$HbA1c < 0))	# 14
sum((lm_proteins_age_HbA1c$Age < 0)&(lm_proteins_age_HbA1c$HbA1c > 0))	# 28
sum((lm_proteins_age_HbA1c$Age < 0)&(lm_proteins_age_HbA1c$HbA1c < 0))	# 2
sum((lm_proteins_age_HbA1c$Age > 0)&(lm_proteins_age_HbA1c$HbA1c > 0)&(lm_proteins_age_HbA1c$BMISignifQ == " SIGNIF"))	# 98
sum((lm_proteins_age_HbA1c$Age > 0)&(lm_proteins_age_HbA1c$HbA1c < 0)&(lm_proteins_age_HbA1c$BMISignifQ == " SIGNIF"))	# 13
sum((lm_proteins_age_HbA1c$Age < 0)&(lm_proteins_age_HbA1c$HbA1c > 0)&(lm_proteins_age_HbA1c$BMISignifQ == " SIGNIF"))	# 27
sum((lm_proteins_age_HbA1c$Age < 0)&(lm_proteins_age_HbA1c$HbA1c < 0)&(lm_proteins_age_HbA1c$BMISignifQ == " SIGNIF"))	# 2



##	Metabolites
ggplot(lm_metab_merged_age_BMI,aes(x=as.numeric(Age),y=as.numeric(BMI),col=HbA1cSignifQ, label=CommonName)) + geom_point() +
  geom_vline(xintercept=c(0,0), col="black") +
  geom_hline(yintercept=c(0,0), col="black") +
  scale_color_manual(values=c("gray70", "palegreen4")) +
  guides(col = FALSE) +
  geom_label_repel(data = subset(lm_metab_merged_age_BMI, subset=(abs(Age)>0.005 & abs(BMI)>0.01 & HbA1c_Q < 0.05)),
                   size = 2, max.overlaps = 7,
                   box.padding   = 0.01, 
                   point.padding = 0.01,
                   segment.color = 'grey50') +
  labs(title="Age and BMI Associations: Metabolite Levels", x ="Age Beta", y = "BMI Beta") +
  theme(title = element_text( size = 14), panel.grid = element_blank(),
        axis.text = element_text( size = 12),
        axis.title = element_text( size = 14))

sum((lm_metab_merged_age_BMI$Age > 0)&(lm_metab_merged_age_BMI$BMI > 0))	# 78
sum((lm_metab_merged_age_BMI$Age > 0)&(lm_metab_merged_age_BMI$BMI < 0))	# 86
sum((lm_metab_merged_age_BMI$Age < 0)&(lm_metab_merged_age_BMI$BMI > 0))	# 39
sum((lm_metab_merged_age_BMI$Age < 0)&(lm_metab_merged_age_BMI$BMI < 0))	# 40
sum((lm_metab_merged_age_BMI$Age > 0)&(lm_metab_merged_age_BMI$BMI > 0)&(lm_metab_merged_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 10
sum((lm_metab_merged_age_BMI$Age > 0)&(lm_metab_merged_age_BMI$BMI < 0)&(lm_metab_merged_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 40
sum((lm_metab_merged_age_BMI$Age < 0)&(lm_metab_merged_age_BMI$BMI > 0)&(lm_metab_merged_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 12
sum((lm_metab_merged_age_BMI$Age < 0)&(lm_metab_merged_age_BMI$BMI < 0)&(lm_metab_merged_age_BMI$HbA1cSignifQ == " SIGNIF"))	# 15


ggplot(lm_metab_merged_age_HbA1c,aes(x=as.numeric(Age),y=as.numeric(HbA1c),col=BMISignifQ, label=CommonName)) + geom_point() +
  geom_vline(xintercept=c(0,0), col="black") +
  geom_hline(yintercept=c(0,0), col="black") +
  scale_color_manual(values=c("gray70", "palegreen4")) +
  guides(col = FALSE) +
  geom_label_repel(data = subset(lm_metab_merged_age_HbA1c, subset=(abs(Age)>0.005 & abs(HbA1c)>0.005 & BMI_Q < 0.01)),
                   size = 2, max.overlaps = 7,
                   box.padding   = 0.01, 
                   point.padding = 0.01,
                   segment.color = 'grey50') +
  labs(title="Age and HbA1c Associations: Metabolite Levels", x ="Age Beta", y = "HbA1c Beta") +
  theme(title = element_text( size = 14), panel.grid = element_blank(),
        axis.text = element_text( size = 12),
        axis.title = element_text( size = 14))

sum((lm_metab_merged_age_HbA1c$Age > 0)&(lm_metab_merged_age_HbA1c$HbA1c > 0))	# 4
sum((lm_metab_merged_age_HbA1c$Age > 0)&(lm_metab_merged_age_HbA1c$HbA1c < 0))	# 49
sum((lm_metab_merged_age_HbA1c$Age < 0)&(lm_metab_merged_age_HbA1c$HbA1c > 0))	# 15
sum((lm_metab_merged_age_HbA1c$Age < 0)&(lm_metab_merged_age_HbA1c$HbA1c < 0))	# 14
sum((lm_metab_merged_age_HbA1c$Age > 0)&(lm_metab_merged_age_HbA1c$HbA1c > 0)&(lm_metab_merged_age_HbA1c$BMISignifQ == " SIGNIF"))	# 4
sum((lm_metab_merged_age_HbA1c$Age > 0)&(lm_metab_merged_age_HbA1c$HbA1c < 0)&(lm_metab_merged_age_HbA1c$BMISignifQ == " SIGNIF"))	# 46
sum((lm_metab_merged_age_HbA1c$Age < 0)&(lm_metab_merged_age_HbA1c$HbA1c > 0)&(lm_metab_merged_age_HbA1c$BMISignifQ == " SIGNIF"))	# 14
sum((lm_metab_merged_age_HbA1c$Age < 0)&(lm_metab_merged_age_HbA1c$HbA1c < 0)&(lm_metab_merged_age_HbA1c$BMISignifQ == " SIGNIF"))	# 13


#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#																				INTERACTION PLOTS
#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
which(anno_genes$GeneName == "CDKN1A")
#[1]	5272
pdf(file = "CDKN1A_Age_Sex.pdf",  width = 8.3, height = 5.8)	
p1<-ggplot(data.frame(Exprs = resid(lm(as.matrix(t(raw_genesE[5272,])) ~ 1 + . , data = covar_Exprsf[c(2:5,7:11)])), Age = covar_Exprsf$Age), aes(x = Age, y = Exprs)) + 
    geom_point(shape=1) + stat_smooth(method = "lm", col = "darkred") + labs(title="CDKN1A", x ="Age", y = "Exprs") + scale_fill_manual(values=c("indianred4","darkred")) +
	theme(aspect.ratio=1, axis.text = element_text(size = 10), axis.title = element_text( size = 12, face = "bold" ), title = element_text( size = 14, face = "italic" ))
p2<-ggplot(data.frame(Exprs = resid(lm(as.matrix(t(raw_genesE[5272,])) ~ 1 + . , data = covar_Exprsf[c(1,3:5,7:11)])), Sex = covar_Exprsf$Sex), aes(x = Sex, y = Exprs, fill=Sex)) + 
    geom_boxplot(alpha=0.9) + labs(x ="Sex", y = "Exprs") + scale_fill_manual(values=c("indianred4","darkred")) + guides(fill="none") +
	theme(aspect.ratio=1, axis.text = element_text(size = 10), axis.title = element_text( size = 12, face = "bold" ), title = element_text( size = 14, face = "italic" ))
p3<-ggplot(data.frame(ExprsP = resid(lm(as.matrix(as.double(raw_genesE[5272,])) ~ 1 + . , data = covar_Exprsf[c(3:5,7:11)])), Age = covar_Exprsf$Age, Sex = covar_Exprsf$Sex), aes(x = Age, y = ExprsP)) + 
	geom_point(shape=1) +
	stat_smooth(method = "lm", col = "darkred") +
	labs(x ="Age", y = "Gene Expression") +
	facet_wrap(~Sex) +
	theme(aspect.ratio=1,
		title = element_text( size = 14, face = "italic" ),
		axis.text = element_text( size = 10),
		axis.title = element_text( size = 12, face = "bold" ),
		strip.text = element_text(size = 12, face = "bold" ))
grid.arrange(p1, p2,p3, layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()

which(anno_proteins$GeneName == "SOST")
#[1]	345
plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[345,4],4],3], sep='')]))
pdf(file = "SOST_Age_Sex.pdf",  width = 8.3, height = 5.8)
plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[4,4],4],3], sep='')]))
p1<-ggplot(data.frame(Exprs = resid(lm(as.matrix(t(raw_proteins_Imp[345,])) ~ 1 + . + plate, data = covar_Protf[c(2:5,7)])), Age = covar_Protf$Age), Age = covar_Protf$Age, aes(x = Age, y = Exprs)) + 
    geom_point(shape=1) + stat_smooth(method = "lm", col = "mediumpurple4") + labs(title="SOST", x ="Age", y = "Protein Levels") + scale_fill_manual(values=c("mediumpurple","mediumpurple4")) +
	theme(aspect.ratio=1, axis.text = element_text(size = 10), axis.title = element_text( size = 12, face = "bold"), title = element_text(size = 14))
p2<-ggplot(data.frame(Exprs = resid(lm(as.matrix(t(raw_proteins_Imp[345,])) ~ 1 + . + plate, data = covar_Protf[c(1,3:5,7)])), Sex = covar_Protf$Sex), aes(x = Sex, y = Exprs, fill=Sex)) + 
    geom_boxplot(alpha=0.9) + labs(x ="Sex", y = "Protein Levels") + scale_fill_manual(values=c("mediumpurple","mediumpurple4")) + guides(fill="none") +
	theme(aspect.ratio=1, axis.text = element_text(size = 10), axis.title = element_text( size = 12, face = "bold" ), title = element_text( size = 14))
p3<-ggplot(data.frame(ExprsP = resid(lm(as.matrix(t(raw_proteins_Imp[345,])) ~ 1 + . + plate, data = covar_Protf[c(3:5,7)])), Age = covar_Protf$Age, Sex = covar_Protf$Sex), aes(x = Age, y = ExprsP)) + 
	geom_point(shape=1) +
	stat_smooth(method = "lm", col = "mediumpurple4") +
	labs(x ="Age", y = "Protein Levels") +
	facet_wrap(~Sex) +
	theme(aspect.ratio=1,
		title = element_text( size = 14, face = "italic" ),
		axis.text = element_text( size = 10),
		axis.title = element_text( size = 12, face = "bold" ),
		strip.text = element_text(size = 12, face = "bold" ))
grid.arrange(p1, p2,p3, layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()

which(anno_met_untarg$BIOCHEMICAL2 == "choline")
#[1]	73
pdf(file = "CholineU_Age_Sex.pdf",  width = 8.3, height = 5.8)
p1<-ggplot(data.frame(Exprs = resid(lm(as.matrix(t(raw_met_untarg_Imp[73,])) ~ 1 + . , data = covar_Untargf[c(2:5,7:8)])), Age = covar_Untargf$Age), aes(x = Age, y = Exprs)) + 
    geom_point(shape=1) + stat_smooth(method = "lm", col = "palegreen4") + labs(title="Choline", x ="Age", y = "Metabolite Levels") + scale_fill_manual(values=c("palegreen3","palegreen4")) +
	theme(aspect.ratio=1, axis.text = element_text(size = 10), axis.title = element_text( size = 12, face = "bold" ), title = element_text( size = 14))
p2<-ggplot(data.frame(Exprs = resid(lm(as.matrix(t(raw_met_untarg_Imp[73,])) ~ 1 + . , data = covar_Untargf[c(1,3:5,7:8)])), Sex = covar_Untargf$Sex), aes(x = Sex, y = Exprs, fill=Sex)) + 
    geom_boxplot(alpha=0.9) + labs(x ="Sex", y = "Metabolite Levels") + scale_fill_manual(values=c("palegreen3","palegreen4")) + guides(fill="none") +
	theme(aspect.ratio=1, axis.text = element_text(size = 10), axis.title = element_text( size = 12, face = "bold" ), title = element_text( size = 14, face = "italic" ))
p3<-ggplot(data.frame(ExprsP = resid(lm(as.matrix(as.double(raw_met_untarg_Imp[73,])) ~ 1 + . , data = covar_Untargf[c(3:5,7:8)])), Age = covar_Untargf$Age, Sex = covar_Untargf$Sex), aes(x = Age, y = ExprsP)) + 
	geom_point(shape=1) +
	stat_smooth(method = "lm", col = "palegreen4") +
	labs(x ="Age", y = "Metabolite Levels") +
	facet_wrap(~Sex) +
	theme(aspect.ratio=1,
		title = element_text( size = 14, face = "italic" ),
		axis.text = element_text( size = 10),
		axis.title = element_text( size = 12, face = "bold" ),
		strip.text = element_text(size = 12, face = "bold" ))
grid.arrange(p1, p2,p3, layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()

which(anno_proteins$GeneName == " LEP")
#[1]	137
plate <- as.factor(as.matrix(covar_Protf[,paste('Plate_',panel[panel[,2] %in% protein_info[paste(protein_info[,3],protein_info[,10],sep="_") %in% anno_proteins[137,4],4],3], sep='')]))
pdf(file = "LEP_BMI_Sex.pdf",  width = 8.3, height = 5.8)	
p1 <- ggplot(data.frame(ExprsP = resid(lm(as.matrix(as.double(raw_proteins_Imp[137,])) ~ 1 + . + plate, data = covar_Protf[c(1,3:5)])), Sex = covar_Protf$Sex), aes(x = Sex, y = ExprsP, fill=Sex)) + 
			geom_boxplot(shape=1) +
			labs(title=" LEP", x ="Sex", y = "Protein Levels") +
			scale_fill_manual(values=c("mediumpurple", "mediumpurple4")) + guides(fill="none") +
			theme(aspect.ratio=1,
				title = element_text( size = 16),
				axis.text = element_text( size = 10),
				axis.title = element_text( size = 16, face = "bold" ),
				strip.text = element_text(size = 12, face = "bold" ))
p2 <- ggplot(data.frame(ExprsP = resid(lm(as.matrix(as.double(raw_proteins_Imp[137,])) ~ 1 + . + plate, data = covar_Protf[c(1:2,4:5)])), BMI = covar_Protf$BMI, Sex = covar_Protf$Sex), aes(x = BMI, y = ExprsP)) + 
			geom_point(shape=1) +
			stat_smooth(method = "lm", col = "mediumpurple4") +
			labs(x ="BMI", y = "Protein Levels") +
			theme(aspect.ratio=1,
				title = element_text( size = 16, face = "italic" ),
				axis.text = element_text( size = 10),
				axis.title = element_text( size = 16, face = "bold" ),
				strip.text = element_text(size = 12, face = "bold" ))
p3 <- ggplot(data.frame(ExprsP = resid(lm(as.matrix(as.double(raw_proteins_Imp[137,])) ~ 1 + . + plate, data = covar_Protf[c(1,4:5)])), BMI = covar_Protf$BMI, Sex = covar_Protf$Sex), aes(x = BMI, y = ExprsP)) + 
			geom_point(shape=1) +
			stat_smooth(method = "lm", col = "mediumpurple4") + 
			labs(x ="BMI", y = "Protein Levels") +
			facet_wrap(~Sex) +
			theme(aspect.ratio=1,
				title = element_text( size = 16, face = "italic" ),
				axis.text = element_text( size = 10),
				axis.title = element_text( size = 16, face = "bold" ),
				strip.text = element_text(size = 12, face = "bold" ))
grid.arrange(p1,p2,p3, layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()

which(anno_met_targ$Metabolite == "Pro")
#	22
pdf(file = "Proline_HbA1c_Sex.pdf",  width = 8.3, height = 5.8)
p1<-ggplot(data.frame(Exprs = resid(lm(as.matrix(t(raw_met_targE[22,])) ~ 1 + . , data = covar_Targf[c(1,3:5,7:9)])), Sex = covar_Targf$Sex), aes(x = Sex, y = Exprs, fill=Sex)) + 
    geom_boxplot(alpha=0.9) + labs(title="Proline", x ="Sex", y = "Metabolite Levels") + scale_fill_manual(values=c("palegreen3","palegreen4")) + guides(fill="none") +
	theme(aspect.ratio=1, axis.text = element_text(size = 10), axis.title = element_text( size = 12, face = "bold" ), title = element_text( size = 14))
p2<-ggplot(data.frame(Exprs = resid(lm(as.matrix(t(raw_met_targE[22,])) ~ 1 + . , data = covar_Targf[c(1:3,5,7:9)])), HbA1c = covar_Targf$HbA1c), aes(x = HbA1c, y = Exprs)) + 
    geom_point(shape=1) + stat_smooth(method = "lm", col = "palegreen4") + labs(x ="HbA1c", y = "Metabolite Levels") + scale_fill_manual(values=c("palegreen3","palegreen4")) +
	theme(aspect.ratio=1, axis.text = element_text(size = 10), axis.title = element_text( size = 12, face = "bold" ), title = element_text( size = 14))
p3<-ggplot(data.frame(ExprsP = resid(lm(as.matrix(as.double(raw_met_targE[22,])) ~ 1 + . , data = covar_Targf[c(1,3,5,7:9)])), HbA1c = covar_Targf$HbA1c, Sex = covar_Targf$Sex), aes(x = HbA1c, y = ExprsP)) + 
	geom_point(shape=1) +
	stat_smooth(method = "lm", col = "palegreen4") +
	labs(x ="HbA1c", y = "Metabolite Levels") +
	facet_wrap(~Sex) +
	theme(aspect.ratio=1,
		title = element_text( size = 14, face = "italic" ),
		axis.text = element_text( size = 10),
		axis.title = element_text( size = 12, face = "bold" ),
		strip.text = element_text(size = 12, face = "bold" ))
grid.arrange(p1, p2,p3, layout_matrix = rbind(c(1, 2), c(3, 3)))
dev.off()