#! /bin/zsh

BASEDIR=$PHYLDATA/rawData/DNASDEnriched/Eukaryortes




for seqfile in `find ./og/*/ -maxdepth 2 -name "*.fas" -type f`; do
    seqfilebasename=`basename $seqfile`
    taxonid=${seqfilebasename%:*}
    taxonname=`grep ^$taxonid$'\t' odb9v1_species.tab | awk '{print $2 "_" $3}'`
    tmp=${seqfilebasename##*:}
    geneid=${tmp%.*}    
    tmp=${seqfile%/*}
    ogname=${tmp##*/}
    header=">"$ogname"|"$taxonname"|"$geneid
    echo $header >> ./og/$ogname.fas
    cat $seqfile >> ./og/$ogname.fas
done

for fastafile in `find ./og/ -maxdepth 1 -name "*.fas" -type f`; do
    fasta_formatter -w 80 -i $fastafile -o .tmp
    mv .tmp $fastafile
done

cd $WRPROJ
