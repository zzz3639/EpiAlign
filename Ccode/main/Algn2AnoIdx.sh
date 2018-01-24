#!/bin/bash

# Usage: Algn2AnoBatch.sh ParaGenomeSearch ParaAnno FileList ExePath
# Files to be processed, should contain "chr{number}a{number}b{number}"
ExePath=$4

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

indexdb=0
while read line ; do
    Names[$indexdb]="$line"
    indexdb=$(($indexdb+1))
done < "${ParaGenome[6]}"

unset chrstr
for ((j=0; j<indexdb; ++j))  
do  
    chrstrline1=${Names[$j]#*_chr}
    chrstrline2=${chrstrline1#*_}
    schrstr=$((${#Names[$j]}-${#chrstrline1}-3))
    lchrstr=$((${#chrstrline1}-${#chrstrline2}+2))
    chrstr[$j]=${Names[$j]:$schrstr:lchrstr}
done

index=0
while read line ; do
    filelist[$index]="$line"
    index=$(($index+1))
done < "$3"

#Processing the files
for ((i=0; i<index; ++i))  
do  
#parsing this filename, find a b chrthis
    strthis=${filelist[$i]}
    strtemp1=${strthis%chr*a*b*}
    spos=${#strtemp1}
    strthis=${strthis:spos}
    strcp=${strthis%%.*}
    strthis="$strcp"
    OLD_IFS=$IFS;
    IFS="ab"
    read -a strtemp2 <<<"${strthis}"
    IFS=$OLD_IFS
    chrstrthis=${strtemp2[0]}
    a=${strtemp2[1]}
    b=${strtemp2[2]}
#find chromosome index of chrstrthis
    for ((j=0; j<indexdb; ++j))  
    do  
        if [ ${chrstrthis} = ${chrstr[$j]} ]
       	then
	    chridx=$j
	    break
	fi
    done
#run the annotation script
    echo 0 >${ParaGenome[10]}ParaAnno.temp
    echo $((${ParaGenome[0]}-1)) >>${ParaGenome[10]}ParaAnno.temp
    echo ${ParaAnno[0]} >>${ParaGenome[10]}ParaAnno.temp
    c=$(($a*${ParaAnno[0]}))
    echo ${chridx} >>${ParaGenome[10]}ParaAnno.temp
    echo $((${ParaAnno[1]}+${c})) >>${ParaGenome[10]}ParaAnno.temp
    echo ${ParaAnno[2]} >>${ParaGenome[10]}ParaAnno.temp
    echo ${ParaGenome[6]} >>${ParaGenome[10]}ParaAnno.temp
    echo >${ParaGenome[10]}Prefix.temp
    ${ExePath}Algn2Ano.out ${ParaGenome[10]}ParaAnno.temp ${ParaGenome[10]}Prefix.temp ${filelist[$i]} >${strtemp1}${chrstr[$chridx]}a${a}b${b}.sseq.anno
done


