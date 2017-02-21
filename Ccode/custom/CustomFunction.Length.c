#include"CustomFunction.h"
#include<math.h>
#define Seg_Min 5

struct CustomOpt
{
    float *num1;
    float *num2;
};

void Custom_Init(unsigned char **sseq1, unsigned short **sseq1_num, float **sseq1_num_align, int *l1, unsigned char **sseq2, unsigned short **sseq2_num, float **sseq2_num_align, int *l2, float *alpha, void **opt, const char *para)
{
    int i;
    (*sseq1_num_align) = (float*)malloc(sizeof(float)*(*l1));
    for(i=0;i<l1;i++){
        if(*(*sseq1_num+i)>Seg_Min)
	    *(*sseq1_num_align+i) = sqrt(1.0/Seg_Min*(*(*sseq1_num+i)));
        else
	    *(*sseq1_num_align+i) = 1.0;
    }
    (*sseq2_num_align) = (float*)malloc(sizeof(float)*(*l2));
    for(i=0;i<l2;i++){
        if(*(*sseq2_num+i)>Seg_Min)
	    *(*sseq2_num_align+i) = sqrt(1.0/Seg_Min*(*(*sseq2_num+i)));
        else
	    *(*sseq2_num_align+i) = 1.0;
    }
    struct CustomOpt *O;
    O = (struct CustomOpt*)malloc(sizeof(struct CustomOpt));
    O->num1 = *sseq1_num_align;
    O->num2 = *sseq2_num_align;
    (*opt) = (void*)O;
    (*alpha) = atof(para);
    return;
}

void Custom_Free(void *opt)
{
    free((struct CustomOpt*)opt);
    return;
}

float Custom_MatchingFunction(unsigned char a, unsigned char b, int i, int j, void *opt)
{
    return (a==b)?2.0:-2.0;
}

float Custom_GapFunction(float alpha, int i, char s, void *opt)
{
    return alpha;
}

