#!/bin/bash

#Usage: ./run.sh Paths_Search_Baseline
index=0
while read line ; do
    MYARRAY[$index]="$line"
    index=$(($index+1))
done < "$1"

index=0
while read line ; do
    Names[$index]="$line"
    index=$(($index+1))
done < "${MYARRAY[0]}"

mkdir ${MYARRAY[6]}

for i in ${Names[@]}; do
    mkdir ${MYARRAY[6]}$i
done


