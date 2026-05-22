#!/usr/bin/bash
#PBS -o ExprsExps.o
#PBS -e ExprsExps.e
#PBS -l nodes=1:ppn=1,mem=20gb,walltime=72:00:00
#PBS -t 1-10643

module load intel/perflibs
module load R/4.0.3
Rscript /users/home/dafmic/CrossSectional/RELATIONS/ExprsExprs.R