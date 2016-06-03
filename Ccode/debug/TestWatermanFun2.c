#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<limits.h>
#include"WatermanFun.h"
#include"StateIO.h"


void readfile(int n, const char *filename, unsigned char *sseq, int *sseq_num)
{
    int i;
    FILE *infile;
    int a,b;
    infile=fopen(filename,"r");
    for(i=0;i<n;i++){
	fscanf(infile,"%d",&a);
	fscanf(infile,"%d",&b);
	sseq[i]=a;
	sseq_num[i]=b;
    }
    fclose(infile);
    return;
}

int main(int argc, char **argv)
{
    int m0;
    int n;
    float alpha;
    int a,b;
    int k;
    FILE *para=fopen(argv[3],"r");
    fscanf(para,"%d",&n);
    fscanf(para,"%d",&m0);
    fscanf(para,"%d",&k);
    fscanf(para,"%f",&alpha);
    fscanf(para,"%d",&a);
    fscanf(para,"%d",&b);
    fclose(para);
    unsigned char *sseq;
    int *sseq_num;
    sseq = (unsigned char *)malloc(sizeof(unsigned char)*n);
    sseq_num = (int *)malloc(sizeof(int)*n);
    readfile(n,argv[1],sseq,sseq_num);

    int m=0;
    int i,j;
    for(i=0;i<n;i++){
	m+=sseq_num[i];
    }

    unsigned char *seq;
    seq = (unsigned char *)malloc(sizeof(unsigned char)*m);
    m=0;
    for(i=0;i<n;i++){
	for(j=0;j<sseq_num[i];j++){
	    seq[m]=sseq[i];
	    m+=1;
	}
    }
    int l;
    l = n;

    struct word_node *bow_seq;

    unsigned char *seq2;
    seq2=(unsigned char *)malloc(sizeof(unsigned char)*(b-a));
    for(i=a;i<b;i++){
	seq2[i-a] = seq[i];
    }
    unsigned char *sseq2;
    unsigned short *sseq2_num;
    int l2;
    l2 = Seq2Sseq(seq2,b-a,&sseq2,&sseq2_num);
    readfile(n,argv[2],sseq,sseq_num);
    m = m0;
    n = m;

    float **map_score;
    unsigned char **map_trace;
    map_score = (float**)malloc(sizeof(float *)*(l+1));
    map_trace = (unsigned char **)malloc(sizeof(unsigned char *)*(l+1));
    for(i=0;i<l+1;i++){
	map_score[i] = (float *)malloc(sizeof(float)*(l2+1));
	map_trace[i] = (unsigned char *)malloc(sizeof(unsigned char)*(l2+1));
    }

    SWA_Even(sseq2, sseq, l2, l, MatchScore_Naive, alpha, NULL, NULL, map_score, map_trace);

    float *maxline;
    maxline = (float *)malloc(sizeof(float)*(l+1));
    for(i=0;i<l+1;i++){
	maxline[i] = 0;
	for(j=0;j<l2+1;j++){
	    maxline[i] = SWF_MAX(map_score[i][j],maxline[i]);
	}
	printf("%f\n",maxline[i]);
    }
    /*
    for(i=0;i<l+m/k;i++){
	maxline[i] = 0;
	for(j=0;j<l2+m/k;j++){
	    maxline[i] = SWF_MAX((map_score[i][j])*(map_trace[i][j]==1),maxline[i]);
	}
    }
    int peaks[3];
    struct pair_node pn[3];
    int p1;
    FILE *alignout;
    int s,t;
    for(i=0;i<3;i++){
	peaks[i] = 0;
	for(j=0;j<l+m/k;j++){
	    if(maxline[peaks[i]]<maxline[j]){
		peaks[i] = j;
	    }
	}
	p1 = 0;
	for(j=0;j<l2+m/k;j++){
            if(map_score[peaks[i]][j]>map_score[peaks[i]][p1])
		p1=j;
	}
	Trace_Bow_Even(map_trace,m/k,p1,peaks[i],pn+i);
	alignout = fopen(argv[3+i],"w");
	Print_Alignment(pn+i,seq2,seq,m,k,alignout);
	fclose(alignout);
        s = pn[i].next->p2;
	t = SWF_MIN(2*peaks[i]-s,l+m/k);
	for(j=s;j<t;j++){
	    maxline[j] = 0.0;
	}
    }
    */
    
    free(seq);
    free(seq2);
    free(sseq);
    free(sseq_num);
    free(maxline);
    return 1;
}






