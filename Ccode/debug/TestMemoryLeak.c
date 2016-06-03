#include<stdio.h>
#include"WatermanFun.h"


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
    int n;
    n = atoi(argv[2]);
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
    n=m;

    struct word_node *bow_seq;
    int l;
    int k=1;
    m = atoi(argv[3]);
    while(1){
	l = Seq2Bow_Slow(seq,n,m,k,&bow_seq);
	for(i=0;i<l;i++){
	    Free_Bow(bow_seq+i);
	}
	free(bow_seq);
    }

    return 1;
}





