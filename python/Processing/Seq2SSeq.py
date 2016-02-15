import sys

state=-1;
statecount=0;
outseq=[];
#read inputseq and process
while True:
    line = sys.stdin.readline();
    if not line:
        break;
    line=str(line);
    print line[0]=='\n';
    if line[0]>='0' and line[0]<='9':
        continue;
    statethis=int(line);
    if statethis==state:
        statecount=statecount+1;
    else:
        if state>=0:
            outseq=outseq+[[state,statecount]];
        statecount=1;
        state=statethis;
#output outseq
for sseqline in outseq:
    sys.stdout.write(str(sseqline[0]));
    sys.stdout.write(' ');
    sys.stdout.write(str(sseqline[1]));
    sys.stdout.write('\n');