#!/usr/bin/bash
#PBS -o BICExprsProt.o
#PBS -e BICExprsProt.e
#PBS -l nodes=1:ppn=1,mem=20gb,walltime=72:00:00
#PBS -t 1-1000

module load intel/perflibs
module load R/4.0.3
Rscript /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_Exprs_Prot.R