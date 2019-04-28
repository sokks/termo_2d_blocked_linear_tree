#!/bin/bash

FILES=data/temp/*
for f in $FILES
do 
    echo "Processing $f file..";
    bin/translate $f $f.txt
done
