#!/bin/bash

#This script sort chromosomes in epigenome indexes(to chr1-22,M,X,Y) for a list of epigenomes
#Usage: batchidxsort.sh paths epigenomeindex ExePath

ExePath=$3

index=0
while read line ; do
    pathsort[$index]="$line"
    index=$(($index+1))
done < "$1"

index=0
while read line ; do
    epigenomeindex[$index]="$line"
    index=$(($index+1))
done < "$2"

for ((i=0; i<index; ++i))  
do  
    pathidxthis=${pathsort[0]}${epigenomeindex[$i]}Files
    ${ExePath}SortChrIndex.out ${pathidxthis} ${pathidxthis}
done 


