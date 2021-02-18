#!/bin/bash
#$ -cwd                         
#$ -j y
#$ -l h_data=8G,h_rt=3:00:00
#$ -o /u/scratch/r/roserao/Outfiles
#$ -e /u/scratch/r/roserao/Outfiles
#!/bin/bash

# update paths
. /u/local/Modules/default/init/modules.sh
module load htslib
module load vcftools
module load plink

# local directories, data directories, results directors
ldir=/u/scratch/r/roserao
ddir=$ldir/data 
rdir=$ldir/data/geno

X=2
while [ $X -le 10 ]
do
    echo $X
    INPUT=$ddir/gene_sample/gene${X}.csv


    tdir=$rdir/$X
        if [ -d $tdir ]; then
            rm -r $tdir
        fi
    mkdir $tdir


    [ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }
    
    OLDIFS=$IFS
    IFS=','


    i=0
    while read gene_id start end chr
    do
        if [ $i != 0 ]; then
            
            #echo "${start}"
            #echo "${end}"
            #echo "${chr}"
            chrm=`echo $chr | sed 's/"//g'`
            gene=`echo $gene_id | sed 's/"//g'`
            echo "##################################"
            echo $gene_id

            gdir=$tdir/$gene
            if [ -d $gdir ]; then
                rm -r $gdir
            fi
            mkdir $gdir
            

            tabix -h http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${chrm}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
                ${chrm}:${start}-${end} > $gdir/${gene}_raw.vcf

            if [ ! -f $gdir/${gene}_raw.vcf ]; then
                echo "extract vcf file failed"
                rm -r $gdir
                continue
            fi

            vcftools --vcf $gdir/${gene}_raw.vcf \
                --keep $ddir/AFR_list.txt \
                --remove-indels --recode --recode-INFO-all \
                --out $gdir/${gene}_AFR
            #rm $gdir/${gene}_AFR.log

            vcftools --vcf  $gdir/${gene}_AFR.recode.vcf --plink \
                --out $gdir/${gene}_AFR.recode
            #rm $gdir/${gene}_AFR.recode.log

            plink --file $gdir/${gene}_AFR.recode --maf 0.05 --make-bed --out $gdir/${gene}_AFR.clean
            rm $gdir/${gene}_AFR.clean.nosex
            #rm $tdir/${gene}_AFR.clean.log

            plink --bfile $gdir/${gene}_AFR.clean --r2 --matrix --out $gdir/${gene}_AFR.clean

            python $ldir/geno_sim.py --gene $gene --size 1000 --dir $gdir
            
        fi

        i=$((i+1))

    done < $INPUT
    
    IFS=$OLDIFS
    X=$((X+1))

done

