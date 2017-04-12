#!/bin/bash

ExePath=$2

index=0
while read line ; do
    ARRAY0[$index]="$line"
    index=$(($index+1))
done < "$1"

mkdir ${ARRAY0[8]}
echo ${ARRAY0[0]} >${ARRAY0[8]}ParaCut.temp
echo ${ARRAY0[1]} >>${ARRAY0[8]}ParaCut.temp
echo ${ARRAY0[6]} >>${ARRAY0[8]}ParaCut.temp
echo ${ARRAY0[7]} >>${ARRAY0[8]}ParaCut.temp
echo ${ARRAY0[8]}Segments/ >>${ARRAY0[8]}ParaCut.temp

${ExePath}CutFolder_Init.sh ${ARRAY0[8]}ParaCut.temp
${ExePath}cut_sseq.out ${ARRAY0[8]}ParaCut.temp

index=0
while read line ; do
    Names[$index]="$line"
    index=$(($index+1))
done < "${ARRAY0[6]}"

for i in ${Names[@]}; do
    mkdir ${ARRAY0[8]}$i
done

echo ${ARRAY0[2]} >${ARRAY0[8]}ParaSearch.temp
echo ${ARRAY0[3]} >>${ARRAY0[8]}ParaSearch.temp
echo ${ARRAY0[4]} >>${ARRAY0[8]}ParaSearch.temp
echo ${ARRAY0[5]} >>${ARRAY0[8]}ParaSearch.temp

for i in ${Names[@]}; do
    echo ${ARRAY0[8]}Segments/${i}.idx >${ARRAY0[8]}PathSearch.temp
    echo ${ARRAY0[8]}Segments/${i}/ >>${ARRAY0[8]}PathSearch.temp
    echo ${ARRAY0[6]} >>${ARRAY0[8]}PathSearch.temp
    echo ${ARRAY0[7]} >>${ARRAY0[8]}PathSearch.temp
    echo ${ARRAY0[8]}${i}/ >>${ARRAY0[8]}PathSearch.temp
    ${ExePath}DatabaseSearchNative.out ${ARRAY0[8]}PathSearch.temp ${ARRAY0[8]}ParaSearch.temp ${ARRAY0[9]} >${ARRAY0[8]}${i}.score
done

