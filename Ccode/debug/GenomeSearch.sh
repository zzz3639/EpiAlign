#!/bin/bash

index=0
while read line ; do
    ARRAY0[$index]="$line"
    index=$(($index+1))
done < "$1"

index=0
while read line ; do
    Names[$index]="$line"
    index=$(($index+1))
done < "${ARRAY0[6]}"

index=0
while read line ; do
    SEG[$index]="$line"
    index=$(($index+1))
done < "${ARRAY0[10]}Segments/${Names[$2]}"


for i in ${SEG[@]}; do
    echo $i
done


