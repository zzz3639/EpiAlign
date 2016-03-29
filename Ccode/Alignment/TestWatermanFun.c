#include<stdio.h>
#include<stdlib.h>
#include"WatermanFun.h"

int main()
{
    unsigned char *seq1;
    unsigned char *seq2;
    int l1=5,l2=10;
    seq1=(unsigned char*)malloc(sizeof(unsigned char)*l1);
    seq2=(unsigned char*)malloc(sizeof(unsigned char)*l2);
    seq1[0]=1; seq1[1]=2; seq1[2]=3; seq1[3]=4; seq1[4]=2;
    seq2[0]=1; seq2[1]=2; seq2[2]=3; seq2[3]=4; seq2[4]=3;
    seq2[5]=2; seq2[6]=1; seq2[7]=2; seq2[8]=3; seq2[9]=4;
    float **map_score;
    unsigned char ** map_trace;
    map_trace=(unsigned char **)malloc(sizeof(unsigned char*)*(l2+1));
    map_score=(float **)malloc(sizeof(float*)*(l2+1));
    int i,j;
    for(i=0;i<l2+1;i++){
	map_trace[i]=(unsigned char *)malloc(sizeof(unsigned char)*(l1+1));
	map_score[i]=(float *)malloc(sizeof(float)*(l1+1));
    }
    SWA_Even(seq1,seq2,5,10,MatchScore_Naive,1.0,NULL,NULL,map_score,map_trace);
    printf("\n\n");
    for(i=0;i<l2;i++){
	for(j=0;j<l1;j++){
	    printf("%d\t",(int)map_score[i+1][j+1]);
	}
	printf("\n");
    }
    return 1;
}


