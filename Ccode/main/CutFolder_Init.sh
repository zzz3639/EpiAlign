#!/bin/bash
index=0
while read line ; do
    MYARRAY[$index]="$line"
    index=$(($index+1))
done < "$1"

index=0
while read line ; do
    Names[$index]="$line"
    index=$(($index+1))
done < "${MYARRAY[2]}"

mkdir ${MYARRAY[4]}

for i in ${Names[@]}; do
    mkdir ${MYARRAY[4]}$i
done


