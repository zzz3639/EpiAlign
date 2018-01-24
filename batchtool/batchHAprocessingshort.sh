#!/bin/bash

#This script automatically processing the epigenome searching result on a list of epigenomes.
#The processing includes annotation to alignment results and concatenation of the alignment scores of all the chromosomes into one file.
#Doing similar thing as batchHAprocessing.sh only top segments are annotated.
#Usage: batchHAprocessingshort.sh batchannoparas Annoparafile epigenomeindex exepath

ExePath=$4

index=0
while read line ; do
    Annoparas[$index]="$line"
    index=$(($index+1))
done < "$1"

index=0
while read line ; do
    epigenomeidx[$index]="$line"
    index=$(($index+1))
done < "$3"

for ((i=0; i<index; ++i))  
do  
    GSparapath=${Annoparas[0]}GS${epigenomeidx[$i]}para

    unset GSparathis
    index1=0
    while read line ; do
        GSparathis[$index1]="$line"
        index1=$(($index1+1))
    done < "${GSparapath}"
    GSoutputpath=${GSparathis[10]}

    unset dbthisindex
    index1=0
    while read line ; do
        dbthisindex[$index1]="$line"
        index1=$(($index1+1))
    done < "${GSparathis[6]}"

    unset chrstr
    for ((j=0; j<index1; ++j))  
    do  
        chrstrline1=${dbthisindex[$j]#*_chr}
        chrstrline2=${chrstrline1#*_}
        schrstr=$(( ${#dbthisindex[$j]}-${#chrstrline1}-3 ))
        lchrstr=$(( ${#chrstrline1}-${#chrstrline2}+2 ))
        chrstr[$j]=${dbthisindex[$j]:$schrstr:lchrstr}
        echo ${chrstr[$j]}
    done
    
    scorefull=${GSoutputpath}${epigenomeidx[$i]}_full.score
    rm ${scorefull}
    touch ${scorefull}
    for ((j=0; j<index1; ++j))  
    do  
        chrthispath=${GSoutputpath}${dbthisindex[$j]}.score
        cat ${chrthispath}|awk -v var=${j} '{print $0"\t"var"\t"NR}' >> ${scorefull}
    done
    
    topkfulltemp=${GSoutputpath}${epigenomeidx[$i]}topkfulltemp
    topknopfulltemp=${GSoutputpath}${epigenomeidx[$i]}topknopfulltemp
    topktemp=${GSoutputpath}${epigenomeidx[$i]}topkidx
    topknoptemp=${GSoutputpath}${epigenomeidx[$i]}topknopidx
    anstemp=${GSoutputpath}anstemp

    lseg=`cat ${scorefull}|wc -l`
    let "wseg=${GSparathis[0]}/${GSparathis[1]}"

    awk '{print $0"\t",($4-$2)/$2}' ${scorefull} >${topkfulltemp}
    sort -k 8 -n -r ${topkfulltemp} >${anstemp}
    head -n ${Annoparas[1]} ${anstemp} | awk '{print $6"\t"$7"\t"$8}' >${topktemp}

    awk '{print $2"\t"$4"\t"$6"\t"$7}' ${scorefull} >${topknopfulltemp}
    ${ExePath}NonOverlapTopK.out ${lseg} ${wseg} ${Annoparas[1]} <${topknopfulltemp} >${topknoptemp}
    topkfolder=${GSoutputpath}${epigenomeidx[$i]}topkfolder
    topknopfolder=${GSoutputpath}${epigenomeidx[$i]}topknopfolder
    mkdir ${topkfolder}
    rm ${topkfolder}/*
    mkdir ${topknopfolder}
    rm ${topknopfolder}/*

    copyfilesidx=${GSoutputpath}${epigenomeidx[$i]}copyfilesidx
    unset copylines
    indexcp=0
    while read line ; do
        copylines[$indexcp]="$line"
        indexcp=$(($indexcp+1))
    done < "${topktemp}"

    rm ${copyfilesidx}
    for ((j=0; j<indexcp; ++j))  
    do  
        read -a copylinethis <<<"${copylines[$j]}"
        copys=$((${GSparathis[1]}*(${copylinethis[1]}-1)))
        copyt=$((${copys}+${GSparathis[0]}))
        copythis=${GSoutputpath}${dbthisindex[${copylinethis[0]}]}/a${copys}b${copyt}
        cp ${copythis}.sseq.algn ${topkfolder}/${chrstr[${copylinethis[0]}]}a${copys}b${copyt}.sseq.algn
        echo ${topkfolder}/${chrstr[${copylinethis[0]}]}a${copys}b${copyt}.sseq.algn >>${copyfilesidx}
    done

    copynopfilesidx=${GSoutputpath}${epigenomeidx[$i]}copynopfilesidx
    unset copylines
    indexcp=0
    while read line ; do
        copylines[$indexcp]="$line"
        indexcp=$(($indexcp+1))
    done < "${topknoptemp}"

    rm ${copynopfilesidx}
    for ((j=0; j<indexcp; ++j))  
    do  
        read -a copylinethis <<<"${copylines[$j]}"
        copys=$((${GSparathis[1]}*(${copylinethis[1]}-1)))
        copyt=$((${copys}+${GSparathis[0]}))
        copythis=${GSoutputpath}${dbthisindex[${copylinethis[0]}]}/a${copys}b${copyt}
        cp ${copythis}.sseq.algn ${topknopfolder}/${chrstr[${copylinethis[0]}]}a${copys}b${copyt}.sseq.algn
        echo ${topknopfolder}/${chrstr[${copylinethis[0]}]}a${copys}b${copyt}.sseq.algn >>${copynopfilesidx}
    done

    ${ExePath}Algn2AnoIdx.sh ${GSparapath} $2 ${copyfilesidx} ${ExePath}
    ${ExePath}Algn2AnoIdx.sh ${GSparapath} $2 ${copynopfilesidx} ${ExePath}

done 



