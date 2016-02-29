import numpy

def motifcount(seq,l,topn):
    n=len(seq);
    strseq=''.join(chr(i) for i in seq);
    #generate motif maps
    mapmotif={};
    for i in range(0,n-l+1):
        keythis=strseq[i:i+l];
        if mapmotif.has_key(keythis):
            mapmotif[keythis]=mapmotif[keythis]+1;
        else:
            mapmotif[keythis]=1;
    #find top motifs
    top=[[0]*topn,['']*topn];
    for keythis in mapmotif.keys():
        if mapmotif[keythis]>top[0][0]:
            top[0][0]=mapmotif[keythis];
            top[1][0]=keythis;
            sortindex=numpy.argsort(top[0]);
            topcount=[0]*topn;
            topmotif=['']*topn;
            for j in range(0, topn):
                topcount[j]=top[0][sortindex[j]];
                topmotif[j]=top[1][sortindex[j]];
            top[0]=topcount;
            top[1]=topmotif;
    #output
    return mapmotif, top;


def randseq(seq):
    n=len(seq);
    #establish the dictionary
    seqsort=sorted(set(seq));
    dic={};
    for i in range(0,len(seqsort)):
        dic[seqsort[i]]=i;
    #compute state frequencies
    fre=[0]*len(seqsort);
    for i in range(0,n):
        fre[dic[seq[i]]]=fre[dic[seq[i]]]+1;
    for i in range(0,len(seqsort)):
        fre[i]=1.0*fre[i]/n;
    #generate randomized sequence
    numpy.random.seed();
    seqrand=[0]*len(seq);
    for i in range(0,n):
        mrand=numpy.random.multinomial(1,fre);
        indexrand=numpy.where(mrand==1);
        indexrand=indexrand[0][0];
        seqrand[i]=seqsort[indexrand];
    #output
    return seqrand;


def randseqgram1(seq):
    n=len(seq);
    #establish the dictionary
    seqsort=sorted(set(seq));
    dic={};
    for i in range(0,len(seqsort)):
        dic[seqsort[i]]=i;
    #compute transition frequency
    fre=[];
    for i in range(0,len(seqsort)):
        fre.append([0]*len(seqsort));
    normfre=[0]*len(seqsort);
    for i in range(1,n):
        fre[dic[seq[i-1]]][dic[seq[i]]]=fre[dic[seq[i-1]]][dic[seq[i]]]+1;
        normfre[dic[seq[i-1]]]=normfre[dic[seq[i-1]]]+1;
    for i in range(0,len(seqsort)):
        for j in range(0,len(seqsort)):
            fre[i][j]=1.0*fre[i][j]/normfre[i];
        normfre[i]=1.0*normfre[i]/(n-1);
    #generate randomized sequence
    numpy.random.seed();
    seqrand=[0]*len(seq);
    mrand=numpy.random.multinomial(1,normfre);
    indexrand=numpy.where(mrand==1);
    indexrand=indexrand[0][0];
    seqrand[0]=seqsort[indexrand];
    for i in range(1,n):
        c1=seqrand[i-1];
        mrand=numpy.random.multinomial(1,fre[dic[c1]]);
        indexrand=numpy.where(mrand==1);
        indexrand=indexrand[0][0];
        seqrand[i]=seqsort[indexrand];
    #output
    return seqrand, fre;



