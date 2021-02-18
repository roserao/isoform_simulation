#!/bin/bash
#$ -cwd                         
#$ -j y
#$ -l h_data=8G,h_rt=6:00:00
#$ -o /u/scratch/r/roserao/Outfiles
#$ -e /u/scratch/r/roserao/Outfiles
#!/bin/bash

# update paths
. /u/local/Modules/default/init/modules.sh
module load python/3.7.2


# local directories, data directories, results directors
ldir=/u/scratch/r/roserao
ddir=$ldir/data
rdir=$ldir/data/geno
sdir=$ldir/data/isoexp_sim

X=2
while [ $X -le 10 ]
do
    echo $X
    tdir=$rdir/$X

    INPUT=$ddir/gene_sample/gene${X}.csv
    [ ! -f $INPUT ] && { echo "$INPUT file not found"; exit 99; }

    OLDIFS=$IFS
    IFS=','

    i=0
    while read gene_id start end chr
    do
        if [ $i != 0 ]; then
            

            chrm=`echo $chr | sed 's/"//g'`
            gene=`echo $gene_id | sed 's/"//g'`
            echo "##################################"
            echo $gene

            gdir=$tdir/$gene
            if [ ! -d $gdir ] 
            then
                echo "Directory /path/to/dir does not exist." 
                continue
            fi

            gsdir=$sdir/$X/$gene
			if [ -d $gsdir ] 
            then
                rm -r $gsdir
            fi
            mkdir $gsdir

            
            for h2g in 0.1 0.2 0.3 0.4 0.5
			do
				gssdir=$sdir/$X/$gene/$h2g
				if [ -d $gssdir ] 
            	then
                	rm -r $gssdir
            	fi
                mkdir $gssdir

  				python3 $ldir/model.py --gene $gene --isoform $X --h2g $h2g --snp_causal 0.1 --snp_ld 2 --sim 100
			done

        fi

        i=$((i+1))

    done < $INPUT
    
    IFS=$OLDIFS
    X=$((X+1))

done

