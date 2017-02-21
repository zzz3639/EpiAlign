#include"CustomFunction.h"
#include"EpiBLASTConstant.h"

#define epsilon 1.5

struct CF_pair
{
    float *p1;
    float *p2;
};

void Custom_Init(unsigned char **sseq1, unsigned short **sseq1_num, float **sseq1_num_align, int *l1, unsigned char **sseq2, unsigned short **sseq2_num, float **sseq2_num_align, int *l2, float *alpha, void **opt, const char *para)
{
    float *f_score;
    f_score = (float*)malloc(sizeof(float)*STATE_MAX);
    int i;
    for(i=0;i<STATE_MAX;i++){
	f_score[i] = 0.0;
    }
    int k;
    k = FILE_CountLine(para,NULL);
    FILE *para_align=fopen(para,"r");
    int t;
    for(i=0;i<k;i++){
	fscanf(para_align,"%d",&t);
	fscanf(para_align,"%f",f_score+t);
    }
    fclose(para_align);
    struct CF_pair *cfp;
    cfp = (struct CF_pair*)malloc(sizeof(struct CF_pair));
    (*alpha) = 1.0;
    (*sseq1_num_align) = (float*)malloc(sizeof(float)*(*l1));
    (*sseq2_num_align) = (float*)malloc(sizeof(float)*(*l2));
    for(i=0;i<*l1;i++){
	*((*sseq1_num_align)+i) = f_score[*((*sseq1)+i)];
    }
    for(i=0;i<*l2;i++){
	*((*sseq2_num_align)+i) = f_score[*((*sseq2)+i)];
    }
    cfp->p1 = *sseq1_num_align;
    cfp->p2 = *sseq2_num_align;
    (*opt) = (void*)cfp;
    free(f_score);
    return;
}

void Custom_Free(void *opt)
{
    free((struct CF_pair*)opt);
    return;
}

float Custom_MatchingFunction(unsigned char a, unsigned char b, int i, int j, void *opt)
{
    return (a==b)?(*(((struct CF_pair*)opt)->p1+i))*2.0:(-(*(((struct CF_pair*)opt)->p1+i))-(*(((struct CF_pair*)opt)->p2+j)))*epsilon;
}

float Custom_GapFunction(float alpha, int i, char s, void *opt)
{
    if(s=='1'){
	return alpha*epsilon*(*((((struct CF_pair*)opt)->p1)+i));
    }
    else{
	return alpha*epsilon*(*((((struct CF_pair*)opt)->p2)+i));
    }
}

