#!/bin/bash

#This script generate parameter files for horizontal alignment by naive comparison method for a list of epigenomes
#Usage: batchHAbaselinepara.sh GSBLparas batchGSBLpaths epigenomeindex exepath

ExePath=$4

indext=0
while read line ; do
    templatein[$indext]="$line"
    indext=$(($indext+1))
done < "$1"

index=0
while read line ; do
    batchBLpaths[$index]="$line"
    index=$(($index+1))
done < "$2"

index=0
while read line ; do
    epigenomeidx[$index]="$line"
    index=$(($index+1))
done < "$3"

for ((i=0; i<index; ++i))  
do  
    parapath=${batchBLpaths[2]}GSBL${epigenomeidx[$i]}para
    epigenomerandpath=${batchBLpaths[0]}${epigenomeidx[$i]}randbatch/
    cat $1 > ${parapath}
    echo ${batchBLpaths[0]}${epigenomeidx[$i]}Files >> ${parapath}
    echo ${batchBLpaths[0]}${epigenomeidx[$i]}/ >> ${parapath}
    echo ${batchBLpaths[0]}${epigenomeidx[$i]}Files >> ${parapath}
    echo ${batchBLpaths[0]}${epigenomeidx[$i]}/ >> ${parapath}
    echo ${epigenomerandpath}${epigenomeidx[$i]}randFiles_1 ${epigenomerandpath}${epigenomeidx[$i]}randFiles_2 ${epigenomerandpath}${epigenomeidx[$i]}randFiles_3 >> ${parapath}
    echo ${epigenomerandpath}${epigenomeidx[$i]}rand_1/ ${epigenomerandpath}${epigenomeidx[$i]}rand_2/ ${epigenomerandpath}${epigenomeidx[$i]}rand_3/ >> ${parapath}
    echo ${batchBLpaths[1]}${epigenomeidx[$i]}output/ >> ${parapath}
done


