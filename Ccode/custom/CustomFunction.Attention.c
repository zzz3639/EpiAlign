#include"CustomFunction.h"
#include"EpiBLASTConstant.h"
#include<stdlib.h>

#define epsilon 1.5

void ComputeScoreFre(unsigned char *sseq, int ns, int w, float *score_fre)
{
    int i,j,k;
    char dicfull[STATE_MAX];
    for(i=0;i<STATE_MAX;i++){
	dicfull[i] = 0;
    }
    for(i=0;i<ns;i++){
	dicfull[sseq[i]] = 1;
    }
    int nc=0;
    for(i=0;i<STATE_MAX;i++){
	if(dicfull[i])
	    nc++;
    }
    int dic[STATE_MAX];
    int *antidic;
    antidic = (int*)malloc(sizeof(int)*nc);
    j=0;
    for(i=0;i<STATE_MAX;i++){
	if(dicfull[i]){
	    antidic[j] = i;
	    dic[i] = j;
	    j++;
	}
    }
    int *idx_temp;
    int *score_temp;
    idx_temp = (int*)malloc(sizeof(int)*nc);
    score_temp = (int*)malloc(sizeof(int)*nc);
    int *sum_temp;
    sum_temp = (int*)malloc(sizeof(int)*nc);
    int t;
    int s1,s2;
    if(ns==0){
    }
    else if(ns<=w+1){
	for(i=0;i<nc;i++){
	    sum_temp[i] = 0;
	}
	for(i=0;i<ns;i++){
	    sum_temp[dic[sseq[i]]]++;
	}
	for(i=0;i<ns;i++){
	    t = dic[sseq[i]];
	    s1=0;
	    s2=0;
	    for(j=0;j<nc;j++){
		if(sum_temp[j]>0){
		    s1++;
		}
		if(sum_temp[j]>sum_temp[t]){
		    s2++;
		}
	    }
	    score_fre[i] = 1.0+1.0*s2/(s1-0.99);
	}
    }
    else if(ns<2*w+1){
	for(i=0;i<nc;i++){
	    sum_temp[i] = 0;
	}
	for(i=0;i<w;i++){
	    sum_temp[dic[sseq[i]]]++;
	}
	for(i=0;i<ns-w;i++){
	    sum_temp[dic[sseq[i+w]]]++;
	    s1=0;
	    s2=0;
	    t = dic[sseq[i]];
	    for(j=0;j<nc;j++){
		if(sum_temp[j]>0){
		    s1++;
		}
		if(sum_temp[j]>sum_temp[t]){
		    s2++;
		}
	    }
	    score_fre[i] = 1.0+1.0*s2/(s1-0.99);
	}
	for(i=ns-w;i<w+1;i++){
	    s1=0;
	    s2=0;
	    t = dic[sseq[i]];
	    for(j=0;j<nc;j++){
		if(sum_temp[j]>0){
		    s1++;
		}
		if(sum_temp[j]>sum_temp[t]){
		    s2++;
		}
	    }
	    score_fre[i] = 1.0+1.0*s2/(s1-0.99);
	}
	for(i=w+1;i<ns;i++){
	    sum_temp[dic[sseq[i-w-1]]]--;
	    s1=0;
	    s2=0;
	    t = dic[sseq[i]];
	    for(j=0;j<nc;j++){
		if(sum_temp[j]>0){
		    s1++;
		}
		if(sum_temp[j]>sum_temp[t]){
		    s2++;
		}
	    }
	    score_fre[i] = 1.0+1.0*s2/(s1-0.99);
	}
    }
    else{
	for(i=0;i<nc;i++){
	    sum_temp[i] = 0;
	}
	for(i=0;i<w;i++){
	    sum_temp[dic[sseq[i]]]++;
	}
        for(i=0;i<w+1;i++){
	    sum_temp[dic[sseq[i+w]]]++;
	    s1=0;
	    s2=0;
	    t = dic[sseq[i]];
	    for(j=0;j<nc;j++){
		if(sum_temp[j]>0){
		    s1++;
		}
		if(sum_temp[j]>sum_temp[t]){
		    s2++;
		}
	    }
	    score_fre[i] = 1.0+1.0*s2/(s1-0.99);
        }
	for(i=w+1;i<ns-w;i++){
	    sum_temp[dic[sseq[i+w]]]++;
	    sum_temp[dic[sseq[i-w-1]]]--;
	    s1=0;
	    s2=0;
	    t = dic[sseq[i]];
	    for(j=0;j<nc;j++){
		if(sum_temp[j]>0){
		    s1++;
		}
		if(sum_temp[j]>sum_temp[t]){
		    s2++;
		}
	    }
	    score_fre[i] = 1.0+1.0*s2/(s1-0.99);
	}
	for(i=ns-w;i<ns;i++){
	    sum_temp[dic[sseq[i-w-1]]]--;
	    s1=0;
	    s2=0;
	    t = dic[sseq[i]];
	    for(j=0;j<nc;j++){
		if(sum_temp[j]>0){
		    s1++;
		}
		if(sum_temp[j]>sum_temp[t]){
		    s2++;
		}
	    }
	    score_fre[i] = 1.0+1.0*s2/(s1-0.99);
	}
    }
    free(sum_temp);
    free(idx_temp);
    free(score_temp);
    free(antidic);
    return;
}

void ComputeScoreGram1(unsigned char *sseq, int ns, int w, int add_score, int *score_sseq)
{
    int i,j;
    int **Link;
    Link = (int**)malloc(sizeof(int*)*STATE_MAX);
    for(i=0;i<STATE_MAX;i++){
	Link[i] = (int*)malloc(sizeof(int)*STATE_MAX);
    }
    for(i=0;i<STATE_MAX;i++){
	for(j=0;j<STATE_MAX;j++){
	    Link[i][j] = -1;
	}
    }

    int *link_sseq_minus;
    link_sseq_minus = (int*)malloc(sizeof(int)*(ns-1));
    for(i=0;i<ns-1;i++){
	link_sseq_minus[i] = Link[sseq[i]][sseq[i+1]];
	Link[sseq[i]][sseq[i+1]] = i;
    }
    int *link_sseq_plus;
    link_sseq_plus = (int*)malloc(sizeof(int)*(ns-1));
    for(i=0;i<ns-1;i++){
	link_sseq_plus[i] = ns-1;
    }
    for(i=0;i<ns-1;i++){
	if(link_sseq_minus[i]<0)
	    continue;
	link_sseq_plus[link_sseq_minus[i]] = i;
    }

    if(ns==0){
    }
    else if(ns==1){
	score_sseq[0] = 1;
    }
    else{
	if(link_sseq_plus[0]>w){
	    score_sseq[0] = add_score;
	}
	else{
	    score_sseq[0] = 1;
	}
	if(ns-2-link_sseq_minus[ns-2]>w){
	    score_sseq[ns-1] = add_score;
	}
	else{
	    score_sseq[ns-1] = 1;
	}
	for(i=1;i<ns-1;i++){
	    score_sseq[i] = 0;
	}
	for(i=1;i<ns-1;i++){
	    if(i-1-link_sseq_minus[i-1]>w){
		score_sseq[i] += add_score;
	    }
	    else{
	    }
            if(link_sseq_plus[i]-i>w){
		score_sseq[i] += add_score;
	    }
	    else{
	    }
	}
	for(i=1;i<ns-1;i++){
	    if(score_sseq[i]==0){
		score_sseq[i] = 1;
	    }
	}
    }
    free(link_sseq_plus);
    free(link_sseq_minus);
    for(i=0;i<STATE_MAX;i++){
	free(Link[i]);
    }
    free(Link);

    return;
}

void RegionalNormalize(float *score_sseq, int ns, int wn, float *final_score_sseq)
{
    int norm_temp;
    float sum_temp;
    float mean_temp;
    int i,j,k;
    if(ns==0){
    }
    else if(ns<=wn+1){
	for(i=0;i<ns;i++){
	    final_score_sseq[i] = score_sseq[i];
	}
    }
    else if(ns<2*wn+1){
	norm_temp = wn;
	sum_temp = 0;
	for(i=0;i<wn;i++){
	    sum_temp += score_sseq[i];
	}
	for(i=0;i<ns-wn;i++){
	    norm_temp++;
	    sum_temp += score_sseq[i+wn];
	    final_score_sseq[i] = 1.0*score_sseq[i]*norm_temp/sum_temp;
	}
	for(i=ns-wn;i<wn+1;i++){
	    final_score_sseq[i] = 1.0*score_sseq[i]*norm_temp/sum_temp;
	}
        for(i=wn+1;i<ns;i++){
	    norm_temp--;
	    sum_temp -= score_sseq[i-wn-1];
	    final_score_sseq[i] = 1.0*score_sseq[i]*norm_temp/sum_temp;
	}
    }
    else{
	norm_temp = wn;
	sum_temp = 0;
	for(i=0;i<wn;i++){
	    sum_temp += score_sseq[i];
	}
	for(i=0;i<wn+1;i++){
	    norm_temp++;
	    sum_temp += score_sseq[i+wn];
	    final_score_sseq[i] = 1.0*score_sseq[i]*norm_temp/sum_temp;
	}
	for(i=wn+1;i<ns-wn;i++){
	    sum_temp += score_sseq[i+wn];
	    sum_temp -= score_sseq[i-wn-1];
	    final_score_sseq[i] = 1.0*score_sseq[i]*norm_temp/sum_temp;
	}
	for(i=ns-wn;i<ns;i++){
	    norm_temp--;
	    sum_temp -= score_sseq[i-wn-1];
	    final_score_sseq[i] = 1.0*score_sseq[i]*norm_temp/sum_temp;
	}
    }
    sum_temp = 0.0;
    for(i=0;i<ns;i++){
	sum_temp += final_score_sseq[i];
    }
    mean_temp = sum_temp/ns;
    for(i=0;i<ns;i++){
	final_score_sseq[i] = 1.0*final_score_sseq[i]/mean_temp;
    }
    return;
}

void AttentionScore(unsigned char *sseq, unsigned short *sseq_num, int ns, int w, int add_score, int wn, float *final_score_sseq)
{
    int i,j,k;
    int *score_sseq;
    score_sseq = (int*)malloc(sizeof(int)*ns);
    ComputeScoreGram1(sseq,ns,w,add_score,score_sseq);
    float *score_fre;
    score_fre = (float*)malloc(sizeof(float)*ns);
    ComputeScoreFre(sseq,ns,w,score_fre);
    float *score_multi;
    score_multi = (float*)malloc(sizeof(float)*ns);
    for(i=0;i<ns;i++){
	score_multi[i] = score_fre[i]*score_sseq[i];
	if(sseq_num[i]<2){
	    score_multi[i] = score_multi[i]*0.5;
	}
    }
    RegionalNormalize(score_multi,ns,wn,final_score_sseq);
    free(score_sseq);
    free(score_fre);
    free(score_multi);
    return;
}

struct CF_pair
{
    float *p1;
    float *p2;
};

void Custom_Init(unsigned char **sseq1, unsigned short **sseq1_num, float **sseq1_num_align, int *l1, unsigned char **sseq2, unsigned short **sseq2_num, float **sseq2_num_align, int *l2, float *alpha, void **opt, const char *para)
{
    int i;
    int k;
    int w;
    int wn;
    int add_score=2;
    FILE *para_align=fopen(para,"r");
    fscanf(para_align,"%d",&w);
    fscanf(para_align,"%d",&wn);
    fclose(para_align);
    struct CF_pair *cfp;
    cfp = (struct CF_pair*)malloc(sizeof(struct CF_pair));
    (*alpha) = 1.0;
    (*sseq1_num_align) = (float*)malloc(sizeof(float)*(*l1));
    (*sseq2_num_align) = (float*)malloc(sizeof(float)*(*l2));
    AttentionScore(*sseq1,*sseq1_num,*l1,w,add_score,wn,*sseq1_num_align);
    for(i=0;i<*l1;i++){
	*((*sseq1_num_align)+i) *= 2.0;
    }
    AttentionScore(*sseq2,*sseq2_num,*l2,w,add_score,wn,*sseq2_num_align);
    for(i=0;i<*l2;i++){
	*((*sseq2_num_align)+i) *= 2.0;
    }
    cfp->p1 = *sseq1_num_align;
    cfp->p2 = *sseq2_num_align;
    (*opt) = (void*)cfp;
    return;
}

void Custom_Free(void *opt)
{
    free((struct CF_pair*)opt);
    return;
}

float Custom_MatchingFunction(unsigned char a, unsigned char b, int i, int j, void *opt)
{
    return (a==b)?(*(((struct CF_pair*)opt)->p1+i)+(*(((struct CF_pair*)opt)->p2+j))):(-(*(((struct CF_pair*)opt)->p1+i))-(*(((struct CF_pair*)opt)->p2+j)))*epsilon;
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

