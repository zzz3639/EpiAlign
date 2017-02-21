#!/bin/bash

index=0
while read line ; do
    ParaGenome[$index]="$line"
    index=$(($index+1))
done < "$1"

index=0
while read line ; do
    ParaAnno[$index]="$line"
    index=$(($index+1))
done < "$2"

index=0
while read line ; do
    Names[$index]="$line"
    index=$(($index+1))
done < "${ParaGenome[6]}"

chrthis=0
for i in ${Names[@]}; do
    frag=0
    while read line ; do
        QueryNames[$frag]="$line"
        frag=$(($frag+1))
    done < "${ParaGenome[10]}Segments/${i}.idx"
    a=0;
    b=${ParaGenome[0]};
    for ((j=0; j<$frag; ++j))  
    do  
        echo 0 >${ParaGenome[10]}ParaAnno.temp
        echo $((${ParaGenome[0]}-1)) >>${ParaGenome[10]}ParaAnno.temp
        echo ${ParaAnno[0]} >>${ParaGenome[10]}ParaAnno.temp
        c=$(($a*${ParaAnno[0]}))
        echo ${chrthis} >>${ParaGenome[10]}ParaAnno.temp
        echo $((${ParaAnno[1]}+${c})) >>${ParaGenome[10]}ParaAnno.temp
        echo ${ParaAnno[2]} >>${ParaGenome[10]}ParaAnno.temp
        echo ${ParaGenome[6]} >>${ParaGenome[10]}ParaAnno.temp
        echo >${ParaGenome[10]}Prefix.temp
        Algn2Ano.out ${ParaGenome[10]}ParaAnno.temp ${ParaGenome[10]}Prefix.temp ${ParaGenome[10]}${i}/a${a}b${b}.sseq.algn >${ParaGenome[10]}${i}/a${a}b${b}.sseq.anno
        a=$(($a+${ParaGenome[1]}))
        b=$(($b+${ParaGenome[1]}))
    done
    chrthis=$(($chrthis+1))  
done


