#include<stdio.h>
#include<stdlib.h>
#include"StateIO.h"
#include"WatermanFun.h"

int main(int argc, char **argv)
{
    if(argc==1){
	printf("\nUsage: run.out w sseq_chromosome1 sseq_chromosome2 sseq_output\n");
	return 1;
    }
    int ns1;
    int ns2;
    unsigned char *sseq1;
    unsigned char *sseq2;
    unsigned short *sseq_num1;
    unsigned short *sseq_num2;
    ns1 = Sseq_ReadFile(argv[2], &sseq1, &sseq_num1, NULL);
    ns2 = Sseq_ReadFile(argv[3], &sseq2, &sseq_num2, NULL);
    int nl1;
    int nl2;
    unsigned char *seq1;
    unsigned char *seq2;
    int nlmix;
    unsigned char *seqmix;
    nl1 = Sseq2Seq(sseq1,sseq_num1,ns1,&seq1);
    nl2 = Sseq2Seq(sseq2,sseq_num2,ns2,&seq2);
    nlmix = nl1;
    seqmix = (unsigned char*)malloc(sizeof(unsigned char)*nlmix);
    int i,k;
    int w;
    w = atoi(argv[1]);
    k = 1;
    for(i=0;i<nl1;i++){
	if(i%w==0){
	    k = 1-k;
	}
        if(k==0){
	    seqmix[i] = seq1[i];
	}
	else{
	    seqmix[i] = seq2[i];
	}
    }
    int nsmix;
    unsigned char *sseqmix;
    unsigned short *sseq_nummix;
    nsmix = Seq2Sseq(seqmix,nlmix,&sseqmix,&sseq_nummix);
    FILE *fout;
    fout = fopen(argv[4],"w");
    for(i=0;i<nsmix;i++){
	fprintf(fout,"%d %d\n",sseqmix[i],sseq_nummix[i]);
    }
    fclose(fout);
    free(sseq1);
    free(sseq_num1);
    free(seq1);
    free(sseq2);
    free(sseq_num2);
    free(seq2);
    free(sseqmix);
    free(sseq_nummix);
    free(seqmix);
    return 1;
}



