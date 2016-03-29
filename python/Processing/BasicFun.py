
def chr_state(i):
    if i>26:
        co=chr((i-27)+ord('A'));
    else:
        co=chr((i-1)+ord('a'));
    return co;

def read_sseq(filename):
    seq=[];
    lenseq=[];
    f=open(filename,'r');
    while True:
        line=f.readline();
        if not line:
            break;
        line=str(line);
        if line[0]>'9' or line[0]<'0':
            continue;
        linethis=map(int,line.split(' '));
        seq.append(linethis[0]);
        lenseq.append(linethis[1]);
    return seq, lenseq;
