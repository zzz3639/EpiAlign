#include<stdio.h>
#include<stdlib.h>
#include"StateIO.h"
#include"WatermanFun.h"
#include"CustomFunction.h"

int main(int argc, char **argv)
{
    if(argc!=5){
	printf("\nError!\nUsage: ./run.out seqfile1 seqfile2 para windowsize\n");
	return 1;
    }
    float alpha;
    int windowsize = atoi(argv[4]);
/*Read Sequences*/
    unsigned char *sseq1;
    unsigned short *sseq1_num;
    int i,j,k;
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
    for(i=0;i<l2p1;i++){
	map_score[i] = (float*)malloc(sizeof(float)*(l1p1));
	map_trace[i] = (unsigned char*)malloc(sizeof(unsigned char)*(l1p1));
    }
    
    SWA_Even(sseq1, sseq2, l1, l2, Custom_MatchingFunction, Custom_GapFunction, alpha, NULL, opt, map_score, map_trace);
    for(i=1;i<l2p1;i++){
	for(j=1;j<l1p1;j++){
	    if(map_trace[i][j]!=1){
		map_score[i][j] = 0.0;
	    }
	    if(sseq1[j-1]!=sseq2[i-1]){
		map_score[i][j] = 0.0;
	    }
	}
    }
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
    /*Print sampled score*/
    int m1,m2;
    m1 = sseq1_cumsum[l1]/windowsize;
    m2 = sseq2_cumsum[l2]/windowsize;
    float **matrixtemp;
    int ki,kj;
    matrixtemp = (float**)malloc(sizeof(float*)*l2p1);
    for(i=0;i<l2p1;i++){
	matrixtemp[i] = (float*)malloc(sizeof(float)*m1);
    }
    float xtemp;
    int t;
    for(i=0;i<l2p1;i++){
	t = 0;
	kj = 1;
	j = windowsize;
	for(;j<=sseq1_cumsum[l1];j+=windowsize){
	    xtemp = map_score[i][kj];
	    while(sseq1_cumsum[kj-1]<j&&kj<l1p1){
		if(map_score[i][kj]>xtemp){
		    xtemp = map_score[i][kj];
		}
		kj++;
	    }
	    matrixtemp[i][t]=xtemp;
	    t++;
	    kj--;
	    if(sseq1_cumsum[kj]==j){
		kj++;
	    }
	}
    }
    float **matrixoutput;
    matrixoutput = (float**)malloc(sizeof(float*)*m2);
    for(i=0;i<m2;i++){
	matrixoutput[i] = (float*)malloc(sizeof(float)*m1);
    }
    t = 0;
    ki = 1;
    i = windowsize;
    float *xtemp2;
    xtemp2 = (float*)malloc(sizeof(float)*m1);
    for(;i<=sseq2_cumsum[l2];i+=windowsize){
	for(j=0;j<m1;j++){
	    xtemp2[j] = matrixtemp[ki][j];
	}
	while(sseq2_cumsum[ki-1]<i&&ki<l2p1){
	    for(j=0;j<m1;j++){
		if(matrixtemp[ki][j]>xtemp2[j]){
		    xtemp2[j] = matrixtemp[ki][j];
		}
	    }
	    ki++;
	}
	for(j=0;j<m1;j++){
	    matrixoutput[t][j] = xtemp2[j];
	}
	t++;
	ki--;
	if(sseq2_cumsum[ki]==i){
	    ki++;
	}
    }
    for(i=0;i<m2;i++){
	for(j=0;j<m1;j++){
	    fprintf(stdout,"%f\t",matrixoutput[i][j]);
	}
        fprintf(stdout,"\n");
    }

/*Free memory*/
    Custom_Free(opt);
    free(xtemp2);
    for(i=0;i<m2;i++){
	free(matrixoutput[i]);
    }
    free(matrixoutput);
    for(i=0;i<l2p1;i++){
	free(matrixtemp[i]);
    }
    free(matrixtemp);
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


