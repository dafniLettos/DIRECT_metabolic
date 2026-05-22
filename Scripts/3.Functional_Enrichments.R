# THIS FILE CONTAINS THE R CODE USED TO PERFORM FUNCTIONAL ENRICHMENT ANALYSES
# ON GENES, PROTEINS AND METABOLITES (TARGETED & UNTARGETED DATASET)
# ASSOCIATED WITH: (A).AGE	(B).SEX	(C).BMI	(D).HBA1C
# FOR GENES WE USED GSEA RANKING BASED ON T-VALUE. ENRICHED GO BIOLOGICAL PROCESSES TERMS WERE USED IN REVIGO TO CREATE A TREEMAP.
# FOR PROTEINS WE USED STRING.db TO PERFORM ENRICHMENT ANALYSES FOR GO BP. ENRICHED TERMS WERE USED IN REVIGO TO CREATE A TREEMAP.
# METABOLITES DATASETS WERE MERGED PRIOR TO ANALYSES BASED ON KEGG ID. ENRISHMENT ANALYSIS WERE PERFORMED USING THE FELLA PACKAGE.

#-------------------------------------------------------------------	LIBRARIES/PACKAGES USED	-------------------------------------------------------------------------
library(tibble)
library(dplyr)
library(ggplot2)
library(rcartocolor)
library(clusterProfiler)
library(pathview)
#library(enrichplot)
library(FELLA)

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------						1. GENE EXPRESSION					---------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lm_genes <- read.delim("lmResultsAdjGenes.txt", header = TRUE)
lm_genes$GeneID2 <- sapply(strsplit(lm_genes$GeneID,"\\."), `[`, 1)
lm_genes <- lm_genes[,c(1:4,41,5:40)]
lm_genes_age <- subset(lm_genes, subset=(Age_Q < 0.05))
lm_genes_sex <- subset(lm_genes, subset=(Sex_Q < 0.05))
lm_genes_bmi <- subset(lm_genes, subset=(BMI_Q < 0.05))
lm_genes_hba1c<- subset(lm_genes, subset=(HbA1c_Q < 0.05))

#		SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)
#		GSE - RANKED ENRICHMENT [t-statistics]:

######################################				AGE			##################################################
gene_list_age_ranked <- lm_genes$Age_Tvalue
names(gene_list_age_ranked) <- lm_genes$GeneID2
#gene_list_age_ranked <- na.omit(gene_list_age_ranked)
gse_genes_age <- gseGO(geneList=sort(gene_list_age_ranked, decreasing = TRUE), 
											ont ="BP", 
											eps = 0,
											keyType = "ENSEMBL", 
											nPermSimple = 10000, 
											minGSSize = 10, 
											maxGSSize = 800, 
											pvalueCutoff = 0.05, 
											verbose = TRUE, 
											OrgDb = "org.Hs.eg.db", 
											pAdjustMethod = "BH")
#require(DOSE)
dotplot(gse_genes_age, showCategory=10, font.size = 12, split=".sign") + facet_grid(.~.sign) + scale_fill_gradient(name = "p.adjust", low = "darkred", high = "indianred") + ggtitle("Enriched GO Biological Processes (GSEA): Age-related Gene Expression")

######################################				SEX			##################################################
gene_list_sex_ranked <- lm_genes$Sex_Tvalue
names(gene_list_sex_ranked) <- lm_genes$GeneID2
gse_genes_sex <- gseGO(geneList=sort(gene_list_sex_ranked, decreasing = TRUE), 
											ont ="BP",
											eps = 0,
											keyType = "ENSEMBL", 
											nPermSimple = 10000, 
											minGSSize = 10, 
											maxGSSize = 800, 
											pvalueCutoff = 0.05, 
											verbose = TRUE, 
											OrgDb = "org.Hs.eg.db", 
											pAdjustMethod = "BH")
dotplot(gse_genes_sex, showCategory=10, font.size = 10, split=".sign") + facet_grid(.~.sign) + scale_fill_gradient(name = "p.adjust", low = "darkred", high = "indianred") + ggtitle("Enriched GO Biological Processes (GSEA): Sex-related Gene Expression")

######################################				BMI			##################################################
gene_list_bmi_ranked <- lm_genes$BMI_Tvalue
names(gene_list_bmi_ranked) <- lm_genes$GeneID2
gse_genes_bmi <- gseGO(geneList=sort(gene_list_bmi_ranked, decreasing = TRUE), 
											ont ="BP",
											eps = 0,
											keyType = "ENSEMBL", 
											nPermSimple = 10000, 
											minGSSize = 10, 
											maxGSSize = 800, 
											pvalueCutoff = 0.05, 
											verbose = TRUE, 
											OrgDb = "org.Hs.eg.db", 
											pAdjustMethod = "BH")
dotplot(gse_genes_bmi, showCategory=10, font.size = 10, split=".sign") + facet_grid(.~.sign) + scale_fill_gradient(name = "p.adjust", low = "darkred", high = "indianred") + ggtitle("Enriched GO Biological Processes (GSEA): BMI-related Gene Expression")

######################################				HbA1c			##################################################
gene_list_hba1c_ranked <- lm_genes$HbA1c_Tvalue
names(gene_list_hba1c_ranked) <- lm_genes$GeneID2
gse_genes_hba1c <- gseGO(geneList=sort(gene_list_hba1c_ranked, decreasing = TRUE), 
											ont ="BP",
											eps = 0,
											keyType = "ENSEMBL", 
											nPermSimple = 10000, 
											minGSSize = 10, 
											maxGSSize = 800, 
											pvalueCutoff = 0.05, 
											verbose = TRUE, 
											OrgDb = "org.Hs.eg.db", 
											pAdjustMethod = "BH")
dotplot(gse_genes_hba1c, showCategory=10, font.size = 10, split=".sign") + facet_grid(.~.sign) + scale_fill_gradient(name = "p.adjust", low = "darkred", high = "indianred") + ggtitle("Enriched GO Biological Processes (GSEA): HbA1c-related Gene Expression")
#----------------------------------------------------------------------------------------------
gse_genes_age_signifGO_GSEA <- as.data.frame(gse_genes_age)
#	1131
gse_genes_sex_signifGO_GSEA <- as.data.frame(gse_genes_sex)
#	466
gse_genes_bmi_signifGO_GSEA <- as.data.frame(gse_genes_bmi)
#	566
gse_genes_hba1c_signifGO_GSEA <- as.data.frame(gse_genes_hba1c)
#	466


#		CREATE TREEMAP USING REVIGO - AGE
##	Revigo Input:
write.table(gse_genes_age_signifGO_GSEA[,c(1,7)], file = "GO_GSEA_Age.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
##	Use "Revigo_BP_TreeMap_Genes_Age.R" to define stuff_ageGenes_freq & create treemap
stuff_ageGenes_freq <- data.frame(1:37,unique(factor(stuff_ageGenes$representative)))
colnames(stuff_ageGenes_freq) <- c("Num","Group")
for(i in 1:nrow(stuff_ageGenes_freq)){
	stuff_ageGenes_freq$Frequency[i] <- sum(stuff_ageGenes$representative == stuff_ageGenes_freq$Group[i])
}
stuff_ageGenes_freq <- stuff_ageGenes_freq[order(stuff_ageGenes_freq[,"Frequency"], decreasing=TRUE),]
stuff_ageGenes_freq$Num <- 1:37
rownames(stuff_ageGenes_freq) <- 1:nrow(stuff_ageGenes_freq)

#	SHOW GROUPS WITH FREQUENCY >1
 ggplot(stuff_ageGenes_freq[1:10,],aes(x=Num, y=Group, size = Frequency, color=Group, label =Group)) +
    geom_point(alpha=0) +
	xlim(-10,20) +
    scale_size(range = c(.1, 13), name="Frequency") +
	guides(color = FALSE, size = FALSE) +
	geom_text(check_overlap = TRUE) +
	scale_color_manual(values = carto_pal(12, "Safe")) +
	theme_void()

#		CREATE TREEMAP USING REVIGO - SEX
##	Revigo Input:
write.table(gse_genes_sex_signifGO_GSEA[,c(1,7)], file = "GO_GSEA_Sex.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
##	Use "Revigo_BP_TreeMap_Genes_BMI.R" to define stuff_sexGenes_freq & create treemap
stuff_sexGenes_freq <- data.frame(1:30,unique(factor(stuff_sexGenes$representative)))
colnames(stuff_sexGenes_freq) <- c("Num","Group")
for(i in 1:nrow(stuff_sexGenes_freq)){
	stuff_sexGenes_freq$Frequency[i] <- sum(stuff_sexGenes$representative == stuff_sexGenes_freq$Group[i])
}
stuff_sexGenes_freq <- stuff_sexGenes_freq[order(stuff_sexGenes_freq[,"Frequency"], decreasing=TRUE),]
stuff_sexGenes_freq$Num <- 1:30
rownames(stuff_sexGenes_freq) <- 1:nrow(stuff_sexGenes_freq)

#	SHOW GROUPS WITH FREQUENCY >1
 ggplot(stuff_sexGenes_freq[1:10,],aes(x=Num, y=Group, size = Frequency, color=Group, label =Group)) +
    geom_point(alpha=0) +
	xlim(-10,20) +
    scale_size(range = c(.1, 17), name="Frequency") +
	guides(color = FALSE, size = FALSE) +
	geom_text(check_overlap = TRUE) +
	scale_color_manual(values = carto_pal(12, "Safe")) +
	theme_void()

##	CREATE TREEMAP USING REVIGO - BMI
##	Revigo Input:
write.table(gse_genes_bmi_signifGO_GSEA[,c(1,7)], file = "GO_GSEA_BMI.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
##	Use "Revigo_BP_TreeMap_Genes_BMI.R" to define stuff_bmiGenes_freq & create treemap
stuff_bmiGenes_freq <- data.frame(1:17,unique(factor(stuff_bmiGenes$representative)))
colnames(stuff_bmiGenes_freq) <- c("Num","Group")
for(i in 1:nrow(stuff_bmiGenes_freq)){
	stuff_bmiGenes_freq$Frequency[i] <- sum(stuff_bmiGenes$representative == stuff_bmiGenes_freq$Group[i])
}
stuff_bmiGenes_freq <- stuff_bmiGenes_freq[order(stuff_bmiGenes_freq[,"Frequency"], decreasing=TRUE),]
stuff_bmiGenes_freq$Num <- 1:17
rownames(stuff_bmiGenes_freq) <- 1:nrow(stuff_bmiGenes_freq)

#	SHOW GROUPS WITH FREQUENCY >1
 ggplot(stuff_bmiGenes_freq[1:10,],aes(x=Num, y=Group, size = Frequency, color=Group, label =Group)) +
    geom_point(alpha=0) +
	xlim(-10,20) +
    scale_size(range = c(.1, 17), name="Frequency") +
	guides(color = FALSE, size = FALSE) +
	geom_text(check_overlap = TRUE) +
	scale_color_manual(values = carto_pal(12, "Safe")) +
	theme_void()

##	CREATE TREEMAP USING REVIGO - HbA1c
##	Revigo Input:
write.table(gse_genes_hba1c_signifGO_GSEA[,c(1,7)], file = "GO_GSEA_HbA1c.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
##	Use "Revigo_BP_TreeMap_Genes_HbA1c.R" to define stuff_hba1cGenes_freq & create treemap
stuff_hba1cGenes_freq <- data.frame(1:29,unique(factor(stuff_hba1cGenes$representative)))
colnames(stuff_hba1cGenes_freq) <- c("Num","Group")
for(i in 1:nrow(stuff_hba1cGenes_freq)){
	stuff_hba1cGenes_freq$Frequency[i] <- sum(stuff_hba1cGenes$representative == stuff_hba1cGenes_freq$Group[i])
}
stuff_hba1cGenes_freq <- stuff_hba1cGenes_freq[order(stuff_hba1cGenes_freq[,"Frequency"], decreasing=TRUE),]
stuff_hba1cGenes_freq$Num <- 1:29
rownames(stuff_hba1cGenes_freq) <- 1:nrow(stuff_hba1cGenes_freq)

#	SHOW GROUPS WITH FREQUENCY >1
 ggplot(stuff_hba1cGenes_freq[1:10,],aes(x=Num, y=Group, size = Frequency, color=Group, label =Group)) +
    geom_point(alpha=0) +
	xlim(-10,20) +
    scale_size(range = c(.1, 17), name="Frequency") +
	guides(color = FALSE, size = FALSE) +
	geom_text(check_overlap = TRUE) +
	scale_color_manual(values = carto_pal(12, "Safe")) +
	theme_void()

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------						2. PROTEINS					-------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lm_proteins <- read.delim("lmResultsAdjProteins.txt", header = TRUE)
lm_proteins$GeneID2 <- trimws(lm_proteins$GeneID2)
lm_proteins$GeneName <- trimws(lm_proteins$GeneName)

######################################				AGE			##################################################

#		ENRICHED TERMS FROM STRING.db ANALYSIS ARE USED AS REVIGO INPUT
##	Use "Revigo_BP_TreeMap_Proteins_Age.R" to define stuff_ageGenes_freq & create treemap
stuff_ageProteins_freq <- data.frame(1:49,unique(factor(stuff_ageProteins$representative)))
colnames(stuff_ageProteins_freq) <- c("Num","Group")
for(i in 1:nrow(stuff_ageProteins_freq)){
	stuff_ageProteins_freq$Frequency[i] <- sum(stuff_ageProteins$representative == stuff_ageProteins_freq$Group[i])
}
stuff_ageProteins_freq <- stuff_ageProteins_freq[order(stuff_ageProteins_freq[,"Frequency"], decreasing=TRUE),]
stuff_ageProteins_freq$Num <- 1:49
rownames(stuff_ageProteins_freq) <- 1:nrow(stuff_ageProteins_freq)

#	SHOW GROUPS WITH HIGH FREQUENCY
 ggplot(stuff_ageProteins_freq[1:12,],aes(x=Num, y=Group, size = Frequency, color=Group, label =Group)) +
    geom_point(alpha=0) +
	xlim(-10,20) +
    scale_size(range = c(.1, 17), name="Frequency") +
	guides(color = FALSE, size = FALSE) +
	geom_text(check_overlap = TRUE) +
	scale_color_manual(values = carto_pal(12, "Safe")) +
	theme_void()
	
######################################				SEX			##################################################

#		ENRICHED TERMS FROM STRING.db ANALYSIS ARE USED AS REVIGO INPUT
##	Use "Revigo_BP_TreeMap_Proteins_Sex.R" to define stuff_sexProteins_freq & create treemap
stuff_sexProteins_freq <- data.frame(1:54,unique(factor(stuff_sexProteins$representative)))
colnames(stuff_sexProteins_freq) <- c("Num","Group")
for(i in 1:nrow(stuff_sexProteins_freq)){
	stuff_sexProteins_freq$Frequency[i] <- sum(stuff_sexProteins$representative == stuff_sexProteins_freq$Group[i])
}
stuff_sexProteins_freq <- stuff_sexProteins_freq[order(stuff_sexProteins_freq[,"Frequency"], decreasing=TRUE),]
stuff_sexProteins_freq$Num <- 1:54
rownames(stuff_sexProteins_freq) <- 1:nrow(stuff_sexProteins_freq)

#	PARENTS PLOT
 ggplot(stuff_sexProteins_freq[1:12,], aes(x=Num, y=Group, size = Frequency, color=Group, label =Group)) +
    geom_point(alpha=0) +
	xlim(-30,30) +
    scale_size(range = c(.1, 10), name="Frequency") +
	guides(color = FALSE, size = FALSE) +
	geom_text(check_overlap = TRUE) +
	scale_color_manual(values = carto_pal(12, "Safe")) +
	theme_void()
	
######################################				BMI			##################################################

#		ENRICHED TERMS FROM STRING.db ANALYSIS ARE USED AS REVIGO INPUT
##	Use "Revigo_BP_TreeMap_Proteins_BMI.R" to define stuff_bmiProteins_freq & create treemap
stuff_bmiProteins_freq <- data.frame(1:33,unique(factor(stuff_bmiProteins$representative)))
colnames(stuff_bmiProteins_freq) <- c("Num","Group")
for(i in 1:nrow(stuff_bmiProteins_freq)){
	stuff_bmiProteins_freq$Frequency[i] <- sum(stuff_bmiProteins$representative == stuff_bmiProteins_freq$Group[i])
}
stuff_bmiProteins_freq <- stuff_bmiProteins_freq[order(stuff_bmiProteins_freq[,"Frequency"], decreasing=TRUE),]
stuff_bmiProteins_freq$Num <- 1:33
rownames(stuff_bmiProteins_freq) <- 1:nrow(stuff_bmiProteins_freq)

#	PARENTS PLOT
 ggplot(stuff_bmiProteins_freq[1:9,], aes(x=Num, y=Group, size = Frequency, color=Group, label =Group)) +
    geom_point(alpha=0) +
	xlim(-30,30) +
    scale_size(range = c(.5, 10), name="Frequency") +
	guides(color = FALSE, size = FALSE) +
	geom_text(check_overlap = TRUE) +
	scale_color_manual(values = carto_pal(12, "Safe")) +
	theme_void()	
	
######################################				HbA1c			##################################################

#		ENRICHED TERMS FROM STRING.db ANALYSIS ARE USED AS REVIGO INPUT
##	Use "Revigo_BP_TreeMap_Proteins_HbA1c.R" to define stuff_hba1cProteins_freq & create treemap
stuff_hba1cProteins_freq <- data.frame(1:30,unique(factor(stuff_hba1cProteins$representative)))
colnames(stuff_hba1cProteins_freq) <- c("Num","Group")
for(i in 1:nrow(stuff_hba1cProteins_freq)){
	stuff_hba1cProteins_freq$Frequency[i] <- sum(stuff_hba1cProteins$representative == stuff_hba1cProteins_freq$Group[i])
}
stuff_hba1cProteins_freq <- stuff_hba1cProteins_freq[order(stuff_hba1cProteins_freq[,"Frequency"], decreasing=TRUE),]
stuff_hba1cProteins_freq$Num <- 1:30
rownames(stuff_hba1cProteins_freq) <- 1:nrow(stuff_hba1cProteins_freq)

#	PARENTS PLOT
 ggplot(stuff_hba1cProteins_freq[1:9,], aes(x=Num, y=Group, size = Frequency, color=Group, label =Group)) +
    geom_point(alpha=0) +
	xlim(-30,30) +
    scale_size(range = c(.5, 10), name="Frequency") +
	guides(color = FALSE, size = FALSE) +
	geom_text(check_overlap = TRUE) +
	scale_color_manual(values = carto_pal(12, "Safe")) +
	theme_void()		

#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------	
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------						3. METABOLITES					----------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#		MERGE METABOLITE DATASETS BASED ON KEGG ID
f <- file.choose()
anno_metab <- read.csv(f)
#	"00_Metabolite_Annotations_TargetedMetabolomics_p150-180.csv"
anno_metabU <- read.delim("Chemical_Annotation_Metabolone.txt", header = TRUE)
lm_metabT <- read.delim("lmResultsAdjTarMet.txt", header = TRUE)
lm_metabU <- read.delim("lmResultsAdjUntarMet.txt", header = TRUE)

#	Add KEGG ID to Targeted Metabolites
for(i in 1:nrow(lm_metabT)){
	if(with(anno_metab, KEGG_ID[R_names == lm_metabT$Metabolite[i]]) != "NA"){
		lm_metabT$KEGG_ID[i] <- with(anno_metab, KEGG_ID[R_names == lm_metabT$Metabolite[i]])}
	else{
		lm_metabT$KEGG_ID[i] <- "NA"
	}
	 print(i)
}
lm_metabT <- lm_metabT[,c(1:3,39,4:38)]

lm_metabT_age <- subset(lm_metabT, subset=(Age_Q < 0.05))
Age_KEGG_IDs_Tar <- as.vector(gsub("\n", ",", lm_metabT_age$KEGG_ID))
Age_KEGG_IDs_Tar <- Age_KEGG_IDs_Tar[nzchar(Age_KEGG_IDs_Tar)]
Age_KEGG_IDs_Tar <- unlist(strsplit(Age_KEGG_IDs_Tar,","))
Age_KEGG_IDs_Tar <- unique(Age_KEGG_IDs_Tar)		## 28 KEGG IDs

##	Add KEGG ID to Untargeted Metabolites
for(i in 1:nrow(lm_metabU)){
	if((with(anno_metabU, KEGG[COMP_ID == lm_metabU$COMP_ID[i]]) != "NA")&&(lm_metabU$COMP_ID[i] %in% anno_metabU$COMP_ID)){
		lm_metabU$KEGG_ID[i] <- with(anno_metabU, KEGG[COMP_ID == lm_metabU$COMP_ID[i]])}
	else{
		lm_metabU$KEGG_ID[i] <- "NA"
	}
	 print(i)
}
lm_metabU <- lm_metabU[,c(1:4,41,5:40)]
lm_metabU$KEGG_ID <- trimws(lm_metabU$KEGG_ID)

######################################				AGE			##################################################

lm_metabU_age <- subset(lm_metabU, subset=(Age_Q < 0.05))
Age_KEGG_IDs_Untar <- as.vector(lm_metabU_age$KEGG_ID)
Age_KEGG_IDs_Untar <- unique(Age_KEGG_IDs_Untar)			## 72-2=70 KEGG IDs

# Create Vector with ALL KEGG IDs (for background)
KEGG_IDs_ALL <- as.vector(gsub("\n", ",", anno_metab$KEGG_ID))
KEGG_IDs_ALL <- KEGG_IDs_ALL[nzchar(KEGG_IDs_ALL)]
KEGG_IDs_ALL <- unlist(strsplit(KEGG_IDs_ALL,","))
KEGG_IDs_ALL <- KEGG_IDs_ALL[nzchar(KEGG_IDs_ALL)]
KEGG_IDs_ALL <- c(KEGG_IDs_ALL, as.vector(anno_metabU$KEGG))
KEGG_IDs_ALL <- unlist(strsplit(KEGG_IDs_ALL,","))
KEGG_IDs_ALL <- KEGG_IDs_ALL[nzchar(KEGG_IDs_ALL)]
KEGG_IDs_ALL <- trimws(KEGG_IDs_ALL)
KEGG_IDs_ALL <- unique(KEGG_IDs_ALL)

# Start with the ENRICHMENT Analysis using FELLA package
# buildGraphFromKEGGREST(organism = "hsa", filter.path = NULL)
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
buildDataFromGraph(keggdata.graph = buildGraphFromKEGGREST(organism = "hsa", filter.path = NULL),
									databaseDir = "meta__hsa_Release_108.0_10_31_Oct_23", internalDir = TRUE, matrices = c("hypergeom", "diffusion", "pagerank"), normality = c("diffusion", "pagerank"),
									dampingFactor = 0.85, niter = 1000)
FELLA <- loadKEGGdata(databaseDir = tail(listInternalDatabases(), 1), internalDir = TRUE, loadMatrix = NULL)	
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
metab_input <- c(Age_KEGG_IDs_Tar, Age_KEGG_IDs_Untar)	# 111
metab_input <- unique(metab_input)		# 103
metab_analysis <- defineCompounds(compounds = metab_input, compoundsBackground = KEGG_IDs_ALL, data = FELLA)
getInput(metab_analysis)
#	80 INPUT METABOLITES
# getExcluded(metab_analysis)
metab_analysis2 <- enrich(compounds = metab_input, compoundsBackground = KEGG_IDs_ALL, method = "hypergeom", data = FELLA)
plot(x = metab_analysis2, 
		method = "hypergeom", 
		main = "My first enrichment using the hypergeometric test in FELLA", 
		threshold = 1, 
		data = FELLA)
metab_analysis2_table <- generateResultsTable(object = metab_analysis2, 
																				method = "hypergeom", 
																				threshold = 1, 
																				data = FELLA)
metab_analysis2_table <- metab_analysis2_table[order(metab_analysis2_table[,"p.value"]),]
metab_analysis2_table$KEGG.name <- sapply(strsplit(metab_analysis2_table$KEGG.name," - "), `[`, 1)

ggplot(metab_analysis2_table, aes(CompoundHits, reorder(KEGG.name, CompoundHits), fill = p.value)) +
    geom_bar(stat = "identity") +
	scale_fill_gradient(name = "Pvalue", low = "palegreen4", high = "palegreen") +
	ggtitle("Enriched KEGG Pathways based on Age-related Metabolites") +
	xlab("Count") +
    theme(
        legend.position = "right",
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_blank())

metab_analysis3 <- enrich(compounds = metab_input, method = "diffusion", approx = "normality", data = FELLA)
metab_analysis3 <- runDiffusion(metab_analysis3, approx = "normality", data = FELLA)
plot(x = metab_analysis3, 
    method = "diffusion", 
    main = "My first enrichment using the diffusion analysis in FELLA", 
    threshold = 0.000001, 
    data = FELLA)
metab_analysis3_table <- generateResultsTable(object = metab_analysis3, 
																				method = "diffusion", 
																				threshold = 1, 
																				data = FELLA)
metab_analysis3_table <- metab_analysis3_table[order(metab_analysis3_table[,"p.score"]),]
metab_analysis3_table$KEGG.name <- sapply(strsplit(metab_analysis3_table$KEGG.name," - "), `[`, 1)
write.table(metab_analysis2_table, file = "./8.ENRICHMENTS/Metabolites_FELLA/AGE/Enrich2.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(metab_analysis3_table, file = "./8.ENRICHMENTS/Metabolites_FELLA/AGE/Diffusion2.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

######################################				SEX			##################################################	

lm_metabT_sex <- subset(lm_metabT, subset=(Sex_Q < 0.05))
Sex_KEGG_IDs_Tar <- as.vector(gsub("\n", ",", lm_metabT_sex$KEGG_ID))
Sex_KEGG_IDs_Tar <- Sex_KEGG_IDs_Tar[nzchar(Sex_KEGG_IDs_Tar)]
Sex_KEGG_IDs_Tar <- unlist(strsplit(Sex_KEGG_IDs_Tar,","))
Sex_KEGG_IDs_Tar <- unique(Sex_KEGG_IDs_Tar)		## 43 KEGG IDs
lm_metabU_sex <- subset(lm_metabU, subset=(Sex_Q < 0.05))
Sex_KEGG_IDs_Untar <- as.vector(lm_metabU_sex$KEGG_ID)
Sex_KEGG_IDs_Untar <- unique(Sex_KEGG_IDs_Untar)			## 98-2=96 KEGG IDs

metab_input_Sex <- c(Sex_KEGG_IDs_Tar, Sex_KEGG_IDs_Untar)	# 120
metab_input_Sex <- unique(metab_input_Sex)		# 111
metab_analysis4 <- defineCompounds(compounds = metab_input_Sex, compoundsBackground = KEGG_IDs_ALL, data = FELLA)
getInput(metab_analysis4)
#	85 INPUT METABOLITES
# getExcluded(metab_analysis)
metab_analysis5 <- enrich(compounds = metab_input_Sex, compoundsBackground = KEGG_IDs_ALL, method = "hypergeom", data = FELLA)
plot(x = metab_analysis5, 
		method = "hypergeom", 
		main = "My first enrichment using the hypergeometric test in FELLA", 
		threshold = 1, 
		data = FELLA)
metab_analysis5_table <- generateResultsTable(object = metab_analysis5, 
																				method = "hypergeom", 
																				threshold = 1, 
																				data = FELLA)
metab_analysis5_table <- metab_analysis5_table[order(metab_analysis5_table[,"p.value"]),]
metab_analysis5_table$KEGG.name <- sapply(strsplit(metab_analysis5_table$KEGG.name," - "), `[`, 1)

ggplot(metab_analysis5_table, aes(CompoundHits, reorder(KEGG.name, CompoundHits), fill = p.value)) +
    geom_bar(stat = "identity") +
	scale_fill_gradient(name = "Pvalue", low = "palegreen4", high = "palegreen") +
	ggtitle("Enriched KEGG Pathways based on Sex-related Metabolites") +
	xlab("Count") +
    theme(
        legend.position = "right",
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_blank())

metab_analysis6 <- enrich(compounds = metab_input_Sex, method = "diffusion", approx = "normality", data = FELLA)
metab_analysis6 <- runDiffusion(metab_analysis6, approx = "normality", data = FELLA)
plot(x = metab_analysis6, 
    method = "diffusion", 
    main = "My first enrichment using the diffusion analysis in FELLA", 
    threshold = 0.000001, 
    data = FELLA)
metab_analysis6_table <- generateResultsTable(object = metab_analysis6, 
																				method = "diffusion", 
																				threshold = 1, 
																				data = FELLA)
metab_analysis6_table <- metab_analysis6_table[order(metab_analysis6_table[,"p.score"]),]
metab_analysis6_table$KEGG.name <- sapply(strsplit(metab_analysis6_table$KEGG.name," - "), `[`, 1)
write.table(metab_analysis5_table, file = "./8.ENRICHMENTS/Metabolites_FELLA/SEX/Enrich2.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(metab_analysis6_table, file = "./8.ENRICHMENTS/Metabolites_FELLA/SEX/Diffusion2.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

######################################				BMI			##################################################	

lm_metabT_bmi <- subset(lm_metabT, subset=(BMI_Q < 0.05))
BMI_KEGG_IDs_Tar <- as.vector(gsub("\n", ",", lm_metabT_bmi$KEGG_ID))
BMI_KEGG_IDs_Tar <- BMI_KEGG_IDs_Tar[nzchar(BMI_KEGG_IDs_Tar)]
BMI_KEGG_IDs_Tar <- unlist(strsplit(BMI_KEGG_IDs_Tar,","))
BMI_KEGG_IDs_Tar <- unique(BMI_KEGG_IDs_Tar)		## 41 KEGG IDs
lm_metabU_bmi <- subset(lm_metabU, subset=(BMI_Q < 0.05))
BMI_KEGG_IDs_Untar <- as.vector(lm_metabU_bmi$KEGG_ID)
BMI_KEGG_IDs_Untar <- unique(BMI_KEGG_IDs_Untar)			## 79-2=77 KEGG IDs

metab_input_BMI <- c(BMI_KEGG_IDs_Tar, BMI_KEGG_IDs_Untar)	# 120
metab_input_BMI <- unique(metab_input_BMI)		# 109
metab_analysis7 <- defineCompounds(compounds = metab_input_BMI, compoundsBackground = KEGG_IDs_ALL, data = FELLA)
getInput(metab_analysis7)
#	82 INPUT METABOLITES
# getExcluded(metab_analysis)
metab_analysis8 <- enrich(compounds = metab_input_BMI, compoundsBackground = KEGG_IDs_ALL, method = "hypergeom", data = FELLA)
plot(x = metab_analysis8, 
		method = "hypergeom", 
		main = "My first enrichment using the hypergeometric test in FELLA", 
		threshold = 1, 
		data = FELLA)
metab_analysis8_table <- generateResultsTable(object = metab_analysis8, 
																				method = "hypergeom", 
																				threshold = 1, 
																				data = FELLA)
metab_analysis8_table <- metab_analysis8_table[order(metab_analysis8_table[,"p.value"]),]
metab_analysis8_table$KEGG.name <- sapply(strsplit(metab_analysis8_table$KEGG.name," - "), `[`, 1)

ggplot(metab_analysis8_table, aes(CompoundHits, reorder(KEGG.name, CompoundHits), fill = p.value)) +
    geom_bar(stat = "identity") +
	scale_fill_gradient(name = "Pvalue", low = "palegreen4", high = "palegreen") +
	ggtitle("Enriched KEGG Pathways based on BMI-related Metabolites") +
	xlab("Count") +
    theme(
        legend.position = "right",
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_blank())

metab_analysis9 <- enrich(compounds = metab_input_BMI, method = "diffusion", approx = "normality", data = FELLA)
metab_analysis9 <- runDiffusion(metab_analysis9, approx = "normality", data = FELLA)
plot(x = metab_analysis9, 
    method = "diffusion", 
    main = "My first enrichment using the diffusion analysis in FELLA", 
    threshold = 0.000001, 
    data = FELLA)
metab_analysis9_table <- generateResultsTable(object = metab_analysis9, 
																				method = "diffusion", 
																				threshold = 1, 
																				data = FELLA)
metab_analysis9_table <- metab_analysis9_table[order(metab_analysis9_table[,"p.score"]),]
metab_analysis9_table$KEGG.name <- sapply(strsplit(metab_analysis9_table$KEGG.name," - "), `[`, 1)
write.table(metab_analysis8_table, file = "./8.ENRICHMENTS/Metabolites_FELLA/BMI/Enrich2.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(metab_analysis9_table, file = "./8.ENRICHMENTS/Metabolites_FELLA/BMI/Diffusion2.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)

######################################				HbA1c			##################################################	

lm_metabT_hba1c <- subset(lm_metabT, subset=(HbA1c_Q < 0.05))
HbA1c_KEGG_IDs_Tar <- as.vector(gsub("\n", ",", lm_metabT_hba1c$KEGG_ID))
HbA1c_KEGG_IDs_Tar <- HbA1c_KEGG_IDs_Tar[nzchar(HbA1c_KEGG_IDs_Tar)]
HbA1c_KEGG_IDs_Tar <- unlist(strsplit(HbA1c_KEGG_IDs_Tar,","))
HbA1c_KEGG_IDs_Tar <- unique(HbA1c_KEGG_IDs_Tar)		## 35 KEGG IDs
lm_metabU_hba1c <- subset(lm_metabU, subset=(HbA1c_Q < 0.05)) #Only 1
HbA1c_KEGG_IDs_Untar <- as.vector(lm_metabU_hba1c$KEGG_ID)
HbA1c_KEGG_IDs_Untar <- unique(HbA1c_KEGG_IDs_Untar)			## 98-2=96 KEGG IDs

metab_input_HbA1c <- c(HbA1c_KEGG_IDs_Tar)	# 35
metab_input_HbA1c  <- unique(metab_input_HbA1c)		# 35
metab_analysis13 <- defineCompounds(compounds = metab_input_HbA1c, compoundsBackground = KEGG_IDs_ALL, data = FELLA)
getInput(metab_analysis13)
#	26 INPUT METABOLITES
# getExcluded(metab_analysis)
metab_analysis10 <- enrich(compounds = metab_input_HbA1c, compoundsBackground = KEGG_IDs_ALL, method = "hypergeom", data = FELLA)
plot(x = metab_analysis10, 
		method = "hypergeom", 
		main = "My first enrichment using the hypergeometric test in FELLA", 
		threshold = 1, 
		data = FELLA)
metab_analysis10_table <- generateResultsTable(object = metab_analysis10, 
																				method = "hypergeom", 
																				threshold = 1, 
																				data = FELLA)
metab_analysis10_table <- metab_analysis10_table[order(metab_analysis10_table[,"p.value"]),]
metab_analysis10_table$KEGG.name <- sapply(strsplit(metab_analysis10_table$KEGG.name," - "), `[`, 1)

ggplot(metab_analysis10_table, aes(CompoundHits, reorder(KEGG.name, CompoundHits), fill = p.value)) +
    geom_bar(stat = "identity") +
	scale_fill_gradient(name = "Pvalue", low = "palegreen4", high = "palegreen") +
	ggtitle("Enriched KEGG Pathways based on HbA1c-related Metabolites") +
	xlab("Count") +
    theme(
        legend.position = "right",
        panel.grid = element_blank(),
        axis.title.y = element_blank(),
        strip.text.x = element_text(size = 14, face = "bold"),
        strip.background = element_blank())

metab_analysis11 <- enrich(compounds = metab_input_HbA1c, method = "diffusion", approx = "normality", data = FELLA)
metab_analysis11 <- runDiffusion(metab_analysis11, approx = "normality", data = FELLA)
plot(x = metab_analysis11, 
    method = "diffusion", 
    main = "My first enrichment using the diffusion analysis in FELLA", 
    threshold = 0.000001, 
    data = FELLA)
metab_analysis11_table <- generateResultsTable(object = metab_analysis11, 
																				method = "diffusion", 
																				threshold = 1, 
																				data = FELLA)
metab_analysis11_table <- metab_analysis11_table[order(metab_analysis11_table[,"p.score"]),]
metab_analysis11_table$KEGG.name <- sapply(strsplit(metab_analysis11_table$KEGG.name," - "), `[`, 1)
write.table(metab_analysis10_table, file = "./8.ENRICHMENTS/Metabolites_FELLA/HbA1c/Enrich2.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
write.table(metab_analysis11_table, file = "./8.ENRICHMENTS/Metabolites_FELLA/HbA1c/Diffusion2.txt", quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = TRUE)
