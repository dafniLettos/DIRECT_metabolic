#!/usr/bin/Rscript
index <- as.integer(Sys.getenv("PBS_ARRAYID"))

#	Load Data
exprs_residuals <- read.table("/users/home/dafmic/CrossSectional/RELATIONS/ExprsResiduals_AgeSUBSET.tab",header = TRUE, sep="\t",stringsAsFactors=FALSE, comment.char = "", check.names = FALSE)
AGE <- read.delim("/users/home/dafmic/CrossSectional/RELATIONS/Exprs_Age3027.tab",header = TRUE)
AGE <- as.numeric(unlist(AGE))

exprsA <- as.matrix(as.double(exprs_residuals[index,3:3029]))
relat_ExprsExprs <- matrix(0,ncol=7,nrow=10643)
colnames(relat_ExprsExprs) <- c("GeneID1","GeneName1","GeneID2","GeneName2","Exprs","Exprs_Tvalue","Exprs_pvalue")

for(j in 1:nrow(exprs_residuals)){
		lmx2 <- lm(exprsA ~ as.matrix(as.double(exprs_residuals[j,3:3029])) + AGE)
		coeff <- summary(lmx2)$coefficients
		relat_ExprsExprs[j,] <- c(as.character(exprs_residuals[i,1:2]), as.character(exprs_residuals[j,1:2]), coeff[2,c(1,3:4)])
}

write.table(relat_ExprsExprs, file = paste("/users/home/dafmic/CrossSectional/RELATIONS/results2/ExprsExprs_Assoc", index, sep="_"), quote = FALSE, sep = "\t ", na = "NA", dec = ".", row.names = FALSE, col.names = FALSE)