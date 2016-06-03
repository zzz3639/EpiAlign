#include<stdio.h>
#include<stdlib.h>
#include"StateIO.h"
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
    struct StateIO_Opt opt;
    opt.s = 0;
    MakeKey_Number(&opt);
    n = FILE_CountLine(argv[1],&opt);
    printf("\n%d\n",n);
    unsigned char *sseq;
    unsigned short *sseq_num;
    opt.s = 0;
    opt.n = 0;
    n = Sseq_ReadFile(argv[1],&sseq,&sseq_num,&opt);

    int i,j;
    int m=0;
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

    int n2;
    n2 = atoi(argv[2]);
    unsigned char *sseq2;
    int *sseq_num2;
    sseq2 = (unsigned char *)malloc(sizeof(unsigned char)*n2);
    sseq_num2 = (int *)malloc(sizeof(int)*n2);
    readfile(n2,argv[1],sseq2,sseq_num2);

    int m2=0;
    for(i=0;i<n2;i++){
	m2 += sseq_num2[i];
    }
    unsigned char *seq2;
    seq2 = (unsigned char *)malloc(sizeof(unsigned char)*m2);
    m2=0;
    for(i=0;i<n2;i++){
	for(j=0;j<sseq_num2[i];j++){
	    seq2[m2] = sseq2[i];
	    m2 += 1;
	}
    }
    n2 = m2;

    printf("\n%d %d\n",n,n2);
    char eq = 1;
    for(i=0;i<n;i++){
	if(seq[i]!=seq2[i]){
	    eq=0;
	    break;
	}
    }
    if(eq==0){
	printf("\nNot Equal!\n");
    }
    else{
	printf("\nEqual.\n");
    }
    return 1;
}


