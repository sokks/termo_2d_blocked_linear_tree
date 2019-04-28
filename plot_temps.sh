#!/bin/bash

FILES=data/temp/*.txt
for f in $FILES
do 
    echo "Processing $f file..";
    num=$(echo $f| cut -c 11-16)
    # echo $SUBSTRING

    # bin/translate $f $f.txt
    python3 vis_2d_nonuniform.py $1 $f data/temp/img/$num.png coolwarm temp
done
