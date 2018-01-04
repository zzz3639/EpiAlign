#include<stdio.h>
#include"StateIO.h"
#include"WatermanFun.h"

int main(int argc, char **argv)
{
    if(argc==1){
	printf("\nUsage: run.out sseq_file_to_be_unified output_file_path\n");
	return 1;
    }
    int ns;
    unsigned char *sseq;
    unsigned short *sseq_num;
    ns = Sseq_ReadFile(argv[1], &sseq, &sseq_num, NULL);
    int nl;
    unsigned char *seq;
    nl = Sseq2Seq(sseq,sseq_num,ns,&seq);
    free(sseq);
    free(sseq_num);
    ns = Seq2Sseq(seq,nl,&sseq,&sseq_num);
    int i;
    FILE *fout;
    fout = fopen(argv[2],"w");
    for(i=0;i<ns;i++){
	fprintf(fout,"%d %d\n",sseq[i],sseq_num[i]);
    }
    fclose(fout);
    free(sseq);
    free(sseq_num);
    free(seq);
    return 1;
}








