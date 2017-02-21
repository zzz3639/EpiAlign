#!/bin/bash

index=0
while read line ; do
    NAMES[$index]="$line"
    index=$(($index+1))
done < "$1"

mkdir $4

rm $3
for i in ${NAMES[@]}; do
    echo ${i%txt.Sseq}rand.Sseq >>$3
done

for i in ${NAMES[@]}; do
    FakeChromosomeGenerator.out $2$i $4${i%txt.Sseq}rand.Sseq
done

