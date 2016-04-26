#include<stdio.h>
#include<stdlib.h>
#include<limits.h>
#include"WatermanFun.h"

int main()
{
    unsigned char *seq1;
    unsigned char *seq2;
    int l1=12,l2=22;
    seq1=(unsigned char*)malloc(sizeof(unsigned char)*l1);
    seq2=(unsigned char*)malloc(sizeof(unsigned char)*l2);
    seq1[0]=seq1[1]=1; seq1[2]=seq1[3]=seq1[4]=2; seq1[5]=seq1[6]=seq1[7]=3; seq1[8]=seq1[9]=4; seq1[10]=seq1[11]=2;
    seq2[0]=seq2[1]=4; seq2[2]=seq2[3]=1; seq2[4]=seq2[5]=2; seq2[6]=seq2[7]=3; seq2[8]=seq2[9]=4;
    seq2[10]=seq2[11]=3; seq2[12]=seq2[13]=2; seq2[14]=seq2[15]=1; seq2[16]=seq2[17]=2; seq2[18]=seq2[19]=3; seq2[20]=seq2[21]=4;
    int i,j;
    printf("\nSeq1=:\n");
    for(i=0;i<l1;i++){
	printf("%d ",(int)seq1[i]);
    }
    printf("\n");
    unsigned char *sseq1;
    unsigned char *sseq2;
    unsigned short *len1;
    unsigned short *len2;
    int sl1,sl2;
    sl1=Seq2Sseq(seq1,l1,&sseq1,&len1);
    for(i=0;i<sl1;i++){
	printf("%d ",(int)sseq1[i]);
    }
    printf("\n");
    for(i=0;i<sl1;i++){
	printf("%d ",(int)len1[i]);
    }
    printf("\n");
    sl2=Seq2Sseq(seq2,l2,&sseq2,&len2);
    for(i=0;i<sl2;i++){
	printf("%d ",(int)sseq2[i]);
    }
    printf("\n");
    for(i=0;i<sl2;i++){
	printf("%d ",(int)len2[i]);
    }
    printf("\n");
    float **map_score;
    unsigned char ** map_trace;
    map_trace=(unsigned char **)malloc(sizeof(unsigned char*)*(l2+1));
    map_score=(float **)malloc(sizeof(float*)*(l2+1));
    for(i=0;i<l2+1;i++){
	map_trace[i]=(unsigned char *)malloc(sizeof(unsigned char)*(l1+1));
	map_score[i]=(float *)malloc(sizeof(float)*(l1+1));
    }
    SWA_Even(seq1,seq2,l1,l2,MatchScore_Naive,1.0,NULL,NULL,map_score,map_trace);
    printf("\n\n");
    for(i=0;i<l2;i++){
	for(j=0;j<l1;j++){
	    printf("%d\t",(int)map_score[i+1][j+1]);
	}
	printf("\n");
    }
    /*Test map malloc*/
    float **map_score2;
    unsigned char **map_trace2;
    Malloc_Map_Compact(sseq1, sseq2, len1, len2, sl1, sl2, &map_score2, &map_trace2);
    SWA_Compact_Even(sseq1, sseq2, len1, len2, sl1, sl2, MatchScore_Naive, 1.0, NULL, map_score2, map_trace2);

    int lp1=sl1+1;
    int lp2=sl2+1;
    int n1=0,n2=0;
    for(i=0;i<sl1;i++)
	n1+=len1[i];
    for(i=0;i<sl2;i++)
	n2+=len2[i];
    int np1=n1+1;
    int np2=n2+1;
    /*Initialize*/
    char *M,*M1;
    M=(char*)malloc(sizeof(char)*np2);
    M1=(char*)malloc(sizeof(char)*np1);
    for(i=0;i<np2;i++)
	M[i]=0;
    for(i=0;i<np1;i++)
	M1[i]=0;
    M[0]=1;
    int k=0;
    for(i=0;i<sl2;i++){
	k+=len2[i];
        M[k]=1;
    }
    M1[0]=1;
    k=0;
    for(i=0;i<sl1;i++){
	k+=len1[i];
        M1[k]=1;
    }
    printf("\n\n\n");
    for(i=1;i<=l2;i++){
	if(M[i]){
	    for(j=1;j<=l1;j++){
		printf("%d\t",(int)map_score2[i][j]);
	    }
	}
	else{
	    int t=1;
	    for(j=1;j<=l1;j++){
                if(M1[j]){
		    printf("%d\t",(int)map_score2[i][t]);
		    t++;
		}
		else{
		    printf("-\t");
		}
	    }
	}
	printf("\n");
    }

    return 1;
}


