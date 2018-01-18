#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"StateIO.h"
#include"WatermanFun.h"

/*if fit_model=="normal", then state frequency is counted on original sequence and compared, if fit_model=="compressed", then state frequency is counted after compression*/
float score_baseline(unsigned char *seq1, int n1, unsigned char *seq2, int n2, char *fit_model, int StateMax)
{
    float ans;
    int i,j,k;
    if(strcmp(fit_model,"normal")==0){
	/*count state frequency*/
	float *StateFre1;
	StateFre1 = (float*)malloc(sizeof(float)*StateMax);
	float *StateFre2;
	StateFre2 = (float*)malloc(sizeof(float)*StateMax);
	for(i=0;i<StateMax;i++){
	    StateFre1[i] = 0;
	    StateFre2[i] = 0;
	}
	for(i=0;i<n1;i++){
	    StateFre1[seq1[i]] += 1.0/n1;
	}
	for(i=0;i<n2;i++){
	    StateFre2[seq2[i]] += 1.0/n2;
	}
	/*count Euclidean distance*/
        ans = 0;
	for(i=0;i<StateMax;i++){
	    ans += (StateFre1[i]-StateFre2[i])*(StateFre1[i]-StateFre2[i]);
	}
	ans = -sqrt(ans);
	free(StateFre1);
	free(StateFre2);
    }
    else if(strcmp(fit_model,"compressed")==0){
        /*count state frequency*/
	float *StateFre1;
	StateFre1 = (float*)malloc(sizeof(float)*StateMax);
	float *StateFre2;
	StateFre2 = (float*)malloc(sizeof(float)*StateMax);
	for(i=0;i<StateMax;i++){
	    StateFre1[i] = 0;
	    StateFre2[i] = 0;
	}
	unsigned char *sseq1;
	unsigned char *sseq2;
	unsigned short *sseq1_num;
	unsigned short *sseq2_num;
	int l1;
	int l2;
	l1 = Seq2Sseq(seq1,n1,&sseq1,&sseq1_num);
	l2 = Seq2Sseq(seq2,n2,&sseq2,&sseq2_num);
	for(i=0;i<l1;i++){
	    StateFre1[sseq1[i]] += 1.0/l1;
	}
	for(i=0;i<l2;i++){
	    StateFre2[sseq2[i]] += 1.0/l2;
	}
	/*count Euclidean distance*/
        ans = 0;
	for(i=0;i<StateMax;i++){
	    ans += (StateFre1[i]-StateFre2[i])*(StateFre1[i]-StateFre2[i]);
	}
	ans = -sqrt(ans);
        free(StateFre1);
	free(StateFre2);
	free(sseq1);
	free(sseq1_num);
	free(sseq2);
	free(sseq2_num);
    }
    else{
	ans=0;
    }
    return ans;
}


int main(int argc, char **argv)
{
    if(argc!=4){
	printf("\nError!\nUsage: ./run.out seqfile1 seqfile2 para\n");
        printf("para a string equals to either ""normal"" or ""compressed""\n  Count the state frequency on sequence either compressed or not\n");
	return 1;
    }
    int i,j,k;
/*Read Sequences*/
    unsigned char *sseq1;
    unsigned short *sseq1_num;
    int l1;
    l1 = Sseq_ReadFile(argv[1],&sseq1,&sseq1_num,NULL);
    unsigned char *sseq2;
    unsigned short *sseq2_num;
    int l2;
    l2 = Sseq_ReadFile(argv[2],&sseq2,&sseq2_num,NULL);
/*decompress*/
    unsigned char *seq1;
    unsigned char *seq2;
    int m1;
    int m2;
    m1 = Sseq2Seq(sseq1, sseq1_num, l1, &seq1);
    m2 = Sseq2Seq(sseq2, sseq2_num, l2, &seq2);
    
/*Count StateMax*/
    char StateMaxArray[STATE_MAX];
    for(i=0;i<STATE_MAX;i++){
	StateMaxArray[i] = 0;
    }

    for(i=0;i<l1;i++){
	StateMaxArray[sseq1[i]] = 1;
    }
    for(i=0;i<l2;i++){
	StateMaxArray[sseq2[i]] = 1;
    }

    int StateMax;
    StateMax = 0;
    for(i=0;i<STATE_MAX;i++){
	if(StateMaxArray[i]==1){
	    StateMax = i+1;
	}
    }

    float ans;
    ans = score_baseline(seq1, m1, seq2, m2, argv[3], StateMax);

    printf("Score=%f\n",ans);
/*Free memory*/
    free(sseq1);
    free(sseq1_num);
    free(sseq2);
    free(sseq2_num);
    free(seq1);
    free(seq2);
    return 1;
}


