#include"CustomFunction.h"

void Custom_Init(unsigned char **sseq1, unsigned short **sseq1_num, float **sseq1_num_align, int *l1, unsigned char **sseq2, unsigned short **sseq2_num, float **sseq2_num_align, int *l2, float *alpha, void **opt, const char *para)
{
    (*opt) = NULL;
    (*alpha) = atof(para);
    (*sseq1_num_align) = (float*)malloc(sizeof(float)*(*l1));
    (*sseq2_num_align) = (float*)malloc(sizeof(float)*(*l2));
    return;
}

void Custom_Free(void *opt)
{
    return;
}

float Custom_MatchingFunction(unsigned char a, unsigned char b, int i, int j, void *opt)
{
    return (a==b)?2.0:-2.0;
}

float Custom_GapFunction(float alpha, int i, void *opt)
{
    return alpha;
}

