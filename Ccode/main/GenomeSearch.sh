#!/bin/bash

index=0
while read line ; do
    ARRAY0[$index]="$line"
    index=$(($index+1))
done < "$1"

mkdir ${ARRAY0[10]}
echo ${ARRAY0[0]} >${ARRAY0[10]}ParaCut.temp
echo ${ARRAY0[1]} >>${ARRAY0[10]}ParaCut.temp
echo ${ARRAY0[6]} >>${ARRAY0[10]}ParaCut.temp
echo ${ARRAY0[7]} >>${ARRAY0[10]}ParaCut.temp
echo ${ARRAY0[10]}Segments/ >>${ARRAY0[10]}ParaCut.temp

CutFolder_Init.sh ${ARRAY0[10]}ParaCut.temp
cut_sseq.out ${ARRAY0[10]}ParaCut.temp

index=0
while read line ; do
    Names[$index]="$line"
    index=$(($index+1))
done < "${ARRAY0[6]}"

for i in ${Names[@]}; do
    mkdir ${ARRAY0[10]}$i
done

echo ${ARRAY0[2]} >${ARRAY0[10]}ParaSearch.temp
echo ${ARRAY0[3]} >>${ARRAY0[10]}ParaSearch.temp
echo ${ARRAY0[4]} >>${ARRAY0[10]}ParaSearch.temp
echo ${ARRAY0[5]} >>${ARRAY0[10]}ParaSearch.temp

for i in ${Names[@]}; do
    echo ${ARRAY0[10]}Segments/${i}.idx >${ARRAY0[10]}PathSearch.temp
    echo ${ARRAY0[10]}Segments/${i}/ >>${ARRAY0[10]}PathSearch.temp
    echo ${ARRAY0[6]} >>${ARRAY0[10]}PathSearch.temp
    echo ${ARRAY0[7]} >>${ARRAY0[10]}PathSearch.temp
    echo ${ARRAY0[8]} >>${ARRAY0[10]}PathSearch.temp
    echo ${ARRAY0[9]} >>${ARRAY0[10]}PathSearch.temp
    echo ${ARRAY0[10]}${i}/ >>${ARRAY0[10]}PathSearch.temp
    DatabaseSearch.out ${ARRAY0[10]}PathSearch.temp ${ARRAY0[10]}ParaSearch.temp ${ARRAY0[11]} >${ARRAY0[10]}${i}.score
done

