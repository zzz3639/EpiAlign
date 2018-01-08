#!/bin/bash

#This script generate parameter files to horizontal alignment for a list of epigenomes
#Usage: batchHApara.sh GStemplate batchGSpara epigenomeindex exepath

ExePath=$4

indext=0
while read line ; do
    templatein[$indext]="$line"
    indext=$(($indext+1))
done < "$1"

index=0
while read line ; do
    batchHApara[$index]="$line"
    index=$(($index+1))
done < "$2"

index=0
while read line ; do
    epigenomeidx[$index]="$line"
    index=$(($index+1))
done < "$3"

for ((i=0; i<index; ++i))  
do  
    parapath=${batchHApara[2]}GS${epigenomeidx[$i]}para
    epirandpath=${batchHApara[0]}${epigenomeidx[$i]}randbatch/
    rm ${parapath}
    touch ${parapath}
    echo ${templatein[0]} >> ${parapath}
    echo ${templatein[1]} >> ${parapath}
    echo ${templatein[2]} >> ${parapath}
    echo ${templatein[3]} >> ${parapath}
    echo ${templatein[4]} >> ${parapath}
    echo ${templatein[5]} >> ${parapath}
    echo ${batchHApara[0]}${epigenomeidx[$i]}Files >> ${parapath}
    echo ${batchHApara[0]}${epigenomeidx[$i]}/ >> ${parapath}
    echo ${epirandpath}${epigenomeidx[$i]}randFiles_1 ${epirandpath}${epigenomeidx[$i]}randFiles_2 ${epirandpath}${epigenomeidx[$i]}randFiles_3 >> ${parapath}
    echo ${epirandpath}${epigenomeidx[$i]}rand_1/ ${epirandpath}${epigenomeidx[$i]}rand_2/ ${epirandpath}${epigenomeidx[$i]}rand_3/ >> ${parapath}
    echo ${batchHApara[1]}${epigenomeidx[$i]}output/ >> ${parapath}
    echo ${templatein[11]} >> ${parapath}
done 



