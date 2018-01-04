#!/bin/bash

ExePath=$2

#read parameters
index=0
while read line ; do
    ARRAY0[$index]="$line"
    index=$(($index+1))
done < "$1"

#read indexes
index1=0
while read line ; do
    Names1[$index1]="$line"
    index1=$(($index1+1))
done < "${ARRAY0[1]}"

index2=0
while read line ; do
    Names2[$index2]="$line"
    index2=$(($index2+1))
done < "${ARRAY0[3]}"

#mix the two epigenomes chromosome by chromosome
mkdir ${ARRAY0[6]}
rm ${ARRAY0[5]}
touch ${ARRAY0[5]}
for ((i=0; i<index1; ++i))  
do  
    echo Mix${Names1[$i]} >>${ARRAY0[5]}
    ${ExePath}mix_chromosomes.out ${ARRAY0[0]} ${ARRAY0[2]}${Names1[$i]} ${ARRAY0[4]}${Names2[$i]} ${ARRAY0[6]}Mix${Names1[$i]}
done 



