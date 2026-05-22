#!/usr/bin/bash
#PBS -o BICProtProt.o
#PBS -e BICProtProt.e
#PBS -l nodes=1:ppn=1,mem=20gb,walltime=72:00:00
#PBS -t 21

module load intel/perflibs
module load R/4.0.3
Rscript /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_Prot_Prot2.R