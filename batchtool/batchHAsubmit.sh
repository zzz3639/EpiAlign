#!/bin/bash

#This script automatically generates submit bash script for a list of epigenomes, uses "screen" command to submit jobs on a server
#Usage: batchHAsubmit.sh paths pathsubmitbash epigenomeindex ExePath

ExePath=$4
submitbashfile=$2

index=0
while read line ; do
    pathsubmit[$index]="$line"
    index=$(($index+1))
done < "$1"

index=0
while read line ; do
    epigenomeindex[$index]="$line"
    index=$(($index+1))
done < "$3"

rm ${submitbashfile}
touch ${submitbashfile}
for ((i=0; i<index; ++i))  
do  
    parapaththis=${pathsubmit[0]}GS${epigenomeindex[$i]}para
    echo screen -dmS GS${epigenomeindex[$i]} ${pathsubmit[1]}GenomeSearch_Path.sh ${parapaththis} ${pathsubmit[1]} >> ${submitbashfile}
done 

chmod u+x ${submitbashfile}

