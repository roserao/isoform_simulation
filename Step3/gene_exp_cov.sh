#!/bin/bash
#$ -cwd                         
#$ -j y
#$ -l h_data=8G,h_rt=3:00:00
#$ -o /u/scratch/r/roserao/Outfiles
#$ -e /u/scratch/r/roserao/Outfiles
#!/bin/bash

# update paths
. /u/local/Modules/default/init/modules.sh
module load R/4.0.2

ldir=/u/scratch/r/roserao

R CMD BATCH $ldir/gene_exp_cov.R
