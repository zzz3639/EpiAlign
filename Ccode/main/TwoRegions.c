#include<stdio.h>
#include<stdlib.h>
#include"StateIO.h"
#include"WatermanFun.h"
#include"CustomFunction.h"

int main(int argc, char **argv)
{
    if(argc!=4){
	printf("\nError!\nUsage: ./run.out seqfile1 seqfile2 para\n");
	return 1;
    }
    float alpha;
/*Read Sequences*/
    unsigned char *sseq1;
    unsigned short *sseq1_num;
    int l1;
    l1 = Sseq_ReadFile(argv[1],&sseq1,&sseq1_num,NULL);
    unsigned char *sseq2;
    unsigned short *sseq2_num;
    int l2;
    l2 = Sseq_ReadFile(argv[2],&sseq2,&sseq2_num,NULL);
/*Do Alignment*/
    void *opt;
/*Initialize before doing alignment*/
    float *sseq1_num_align;
    float *sseq2_num_align;
    Custom_Init(&sseq1,&sseq1_num,&sseq1_num_align,&l1,&sseq2,&sseq2_num,&sseq2_num_align,&l2,&alpha,&opt,argv[3]);
    int l1p1,l2p1;
    l1p1 = l1+1; l2p1 = l2+1;
    float **map_score;
    unsigned char **map_trace;
    map_score = (float**)malloc(sizeof(float*)*(l2p1));
    map_trace = (unsigned char**)malloc(sizeof(unsigned char*)*(l2p1));
    int i,j,k;
    for(i=0;i<l2p1;i++){
	map_score[i] = (float*)malloc(sizeof(float)*(l1p1));
	map_trace[i] = (unsigned char*)malloc(sizeof(unsigned char)*(l1p1));
    }
    
    SWA_Even(sseq1, sseq2, l1, l2, Custom_MatchingFunction, Custom_GapFunction, alpha, NULL, opt, map_score, map_trace);
/*Output result to stdout*/
    /*Find maximal score*/
    int *sseq1_cumsum;
    sseq1_cumsum = (int*)malloc(sizeof(int)*(l1+1));
    sseq1_cumsum[0] = 0;
    for(i=0;i<l1;i++){
	sseq1_cumsum[i+1] = sseq1_cumsum[i]+sseq1_num[i];
    }
    int *sseq2_cumsum;
    sseq2_cumsum = (int*)malloc(sizeof(int)*(l2+1));
    sseq2_cumsum[0] = 0;
    for(i=0;i<l2;i++){
	sseq2_cumsum[i+1] = sseq2_cumsum[i]+sseq2_num[i];
    }
    int p1,p2;
    p1=0;p2=0;
    for(i=0;i<l2p1;i++){
	for(j=0;j<l1p1;j++){
	    if(map_score[i][j]>map_score[p2][p1]){
		p1 = j;
		p2 = i;
	    }
	}
    }
    /*Do traceback*/
    struct pair_node align_this;
    Trace_Even(map_trace,p1,p2,&align_this);
    printf("seq1:%d-%d, seq2:%d-%d\n",sseq1_cumsum[align_this.next->p1],sseq1_cumsum[p1]-1,sseq2_cumsum[align_this.next->p2],sseq2_cumsum[p2]-1);
    printf("Score=%f\n",map_score[p2][p1]);
    Print_Alignment_Sseq_Even(&align_this, sseq1, sseq2, sseq1_num, sseq2_num, stdout);


/*Free memory*/
    Custom_Free(opt);
    Free_Alignment(&align_this);
    for(i=0;i<l2p1;i++){
	free(map_score[i]);
	free(map_trace[i]);
    }
    free(map_score);
    free(map_trace);
    free(sseq1);
    free(sseq1_num);
    free(sseq2);
    free(sseq2_num);
    free(sseq1_cumsum);
    free(sseq2_cumsum);
    free(sseq1_num_align);
    free(sseq2_num_align);
    return 1;
}


