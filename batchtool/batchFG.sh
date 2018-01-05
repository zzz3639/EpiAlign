#!/bin/bash

#This script generate fake epigenomes for a list of epigenomes
#Usage: batchFG.sh paths epigenomeindex ExePath

ExePath=$3

index=0
while read line ; do
    pathFG[$index]="$line"
    index=$(($index+1))
done < "$1"

index=0
while read line ; do
    epigenomeindex[$index]="$line"
    index=$(($index+1))
done < "$2"

for ((i=0; i<index; ++i))  
do  
    pathidxthis=${pathFG[0]}${epigenomeindex[$i]}Files
    pathdbthis=${pathFG[0]}${epigenomeindex[$i]}/
    randbatchpath=${pathFG[1]}${epigenomeindex[$i]}randbatch/
    mkdir ${randbatchpath}
    ${ExePath}FakeGenomeGenerator.sh ${pathidxthis} ${pathdbthis} ${randbatchpath}${epigenomeindex[$i]}randFiles_1 ${randbatchpath}${epigenomeindex[$i]}rand_1/ ${ExePath}
    ${ExePath}FakeGenomeGenerator.sh ${pathidxthis} ${pathdbthis} ${randbatchpath}${epigenomeindex[$i]}randFiles_2 ${randbatchpath}${epigenomeindex[$i]}rand_2/ ${ExePath}
    ${ExePath}FakeGenomeGenerator.sh ${pathidxthis} ${pathdbthis} ${randbatchpath}${epigenomeindex[$i]}randFiles_3 ${randbatchpath}${epigenomeindex[$i]}rand_3/ ${ExePath}
done 




