#!/bin/bash

#Run horizontal alignment by naive comparison method
#Usage: ./DatabaseSearch_baseline.sh ParaPaths_Search_Baseline.example ExePath

ExePath=$2

index=0
while read line ; do
    MYARRAY[$index]="$line"
    index=$(($index+1))
done < "$1"

mkdir ${MYARRAY[12]}

Parafile=${MYARRAY[12]}Para_Search_Baseline.temp
Pathsfile=${MYARRAY[12]}Paths_Search_Baseline.temp

head -n 6 $1 >${Parafile}
tail -n +7 $1 >${Pathsfile}

${ExePath}DatabaseSearch_baseline_mkfolder.sh ${Pathsfile}
${ExePath}DatabaseSearch_baseline.out ${Pathsfile} ${Parafile}


