#!/usr/bin/bash
#PBS -o BICUntarMetTarMet.o
#PBS -e BICUntarMetTarMet.e
#PBS -l nodes=1:ppn=1,mem=20gb,walltime=72:00:00
#PBS -t 1-10

module load intel/perflibs
module load R/4.0.3
Rscript /users/home/dafmic/CrossSectional/RELATIONS/BIC/Scripts/BIC_UntarMet_UntarMet.R