#! /bin/zsh

for file in `find $PHYLDATA/modified/ -name "*.afa" -maxdepth 3`; do
    tmp=${file#*modified/}
    dir=${tmp%/*}
    filename=${tmp##*/}
    raxmlHPC -f a -s $file -n $filename -w $PHYLPROJ/phylinference/$dir/ -m GTRCAT -p 12345 -O -x 12345 -# 200
done
