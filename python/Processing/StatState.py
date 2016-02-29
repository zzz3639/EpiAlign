import StatStateFun as SF
import sys
import shelve
#Usage: python StatState.py MotifLength OutFileName <ShortSequenceFile

topn=200;
l=int(sys.argv[1]);
#read seq and len
seq=[];
lenseq=[];
while True:
    line = sys.stdin.readline();
    if not line:
        break;
    line=str(line);
    if line[0]>'9' or line[0]<'0':
        continue;
    linethis=map(int,line.split(' '));
    seq.append(linethis[0]);
    lenseq.append(linethis[1]);
strseq=''.join(chr(i) for i in seq);
#establish the dictionary
seqsort=sorted(set(seq));
quies=seqsort[-1];
dic={};
for i in range(0,len(seqsort)):
    dic[seqsort[i]]=i;
#count average segment length
n=len(seq);
avelen_quies=1.0*reduce(lambda x,y:x+y,lenseq)/n;
avelen_noquies=0.0;
num_noquies=0;
for i in range(0,n):
    if seq[i]!=quies:
        avelen_noquies=avelen_noquies+lenseq[i];
        num_noquies=num_noquies+1;
avelen_noquies=1.0*avelen_noquies/num_noquies;
#count motif frequency
mapmotif_raw, top_raw = SF.motifcount(seq,l,topn);
mapmotif_rand, top_rand = SF.motifcount(SF.randseq(seq),l,topn);
#count baseline
mapmotif_gram1, top_gram1 = SF.motifcount(SF.randseqgram1(seq),l,topn);
#output
outfilename=sys.argv[2];
shelve_save=shelve.open(outfilename,'n');
shelve_save['mapmotif_raw']=mapmotif_raw;
shelve_save['mapmotif_rand']=mapmotif_rand;
shelve_save['mapmotif_gram1']=mapmotif_gram1;
shelve_save['top_raw']=top_raw;
shelve_save['top_rand']=top_rand;
shelve_save['top_gram1']=top_gram1;
shelve_save['avelen_quies']=avelen_quies;
shelve_save['avelen_noquies']=avelen_noquies;
shelve_save.close();