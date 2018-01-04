#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"WatermanFun.h"
#include"StateIO.h"
#include"CustomFunction.h"

/*in case of very long line for randomized baseline*/
#define n_rand_max 5

void Remove_Peak(float **maxline, int *L, int h, int p, int l2, int dist_peak_l, int dist_peak_u)
{
    int l;
    int u;
    l = SWF_MAX(p-l2*dist_peak_l,0);
    u = SWF_MIN(p+l2*dist_peak_u,L[h]);
    int i;
    for(i=l;i<=u;i++){
	maxline[h][i] = 0.0;
    }
    return;
}

void Find_Peak(float **maxline, int n, int *L, int *chr, int *pos)
{
    int h;
    int p;
    h = 0;
    p = 0;
    int i,j;
    for(i=0;i<n;i++){
	for(j=0;j<=L[i];j++){
	    if(maxline[i][j]>maxline[h][p]){
		h = i;
		p = j;
	    }
	}
    }
    *chr = h;
    *pos = p;
    return;
}

int main(int argc, char **argv)
{
    if(argc<2){
	printf("\nUsage:\n./run Paths_Search Para_Search Para_align\n");
	return 1;
    }
/*read searching parameters from file*/
    FILE *para_file=fopen(argv[2],"r");
    int dist_peak_l;
    int dist_peak_u;
    int record_peak;
    float th_peak;
    fscanf(para_file,"%d",&dist_peak_l);
    fscanf(para_file,"%d",&dist_peak_u);
    fscanf(para_file,"%d",&record_peak);
    fscanf(para_file,"%f",&th_peak);
    fclose(para_file);
/*read paths*/
    char Idxq[FileNameLength_MAX];
    char Pq[FileNameLength_MAX];
    char Idxf[FileNameLength_MAX];
    char Pf[FileNameLength_MAX];
    char *Idxr;
    char *Pr;
    Idxr = (char*)malloc(sizeof(char)*FileNameLength_MAX*n_rand_max);
    Pr = (char*)malloc(sizeof(char)*FileNameLength_MAX*n_rand_max);
    char Po[FileNameLength_MAX];
    FILE *seq_paths=fopen(argv[1],"r");
    fscanf(seq_paths,"%s",Idxq);
    fscanf(seq_paths,"%s",Pq);
    fscanf(seq_paths,"%s",Idxf);
    fscanf(seq_paths,"%s",Pf);
    fgetc(seq_paths);
    fgets(Idxr,FileNameLength_MAX*n_rand_max,seq_paths);
    fgets(Pr,FileNameLength_MAX*n_rand_max,seq_paths);
    fscanf(seq_paths,"%s",Po);
    fclose(seq_paths);
    int i,j;
    int n_rand;
    /*count number of randomized baseline datasets*/
    j=1;
    char c;
    i=0;
    while((c=Idxr[i])!='\0'){
	if(c==' '){
	    j++;
	}
	i++;
    }
    n_rand = j;
    /*replace '\n' by '\0'*/
    i=0;
    while(Idxr[i]!='\0'){
        if(Idxr[i]=='\n'){
	    Idxr[i] = '\0';
	    break;
	}
	i++;
    }
    i=0;
    while(Pr[i]!='\0'){
	if(Pr[i]=='\n'){
	    Pr[i] = '\0';
	    break;
	}
	i++;
    }
    /*semitation*/
    char **Idxr_names;
    char **Pr_names;
    Idxr_names = (char**)malloc(sizeof(char*)*j);
    Pr_names = (char**)malloc(sizeof(char*)*j);
    for(i=0;i<n_rand;i++){
	Idxr_names[i] = (char*)malloc(sizeof(char)*FileNameLength_MAX);
	Pr_names[i] = (char*)malloc(sizeof(char)*FileNameLength_MAX);
    }
    int k;
    i=0; j=0; k=0;
    while((c=Idxr[i])!='\0'){
	Idxr_names[j][i-k] = Idxr[i];
	if(c==' '){
	    Idxr_names[j][i-k] = '\0';
	    k = i+1;
	    j++;
	}
	i++;
    }
    i=0; j=0; k=0;
    while((c=Pr[i])!='\0'){
	Pr_names[j][i-k] = Pr[i];
	if(c==' '){
	    Pr_names[j][i-k] = '\0';
	    k = i+1;
	    j++;
	}
	i++;
    }
    

/*read filenames*/
    int nf;
    char **filenames;
    nf = Lines_ReadFile(Idxf,&filenames,NULL);

    int *nr;
    char ***filenamesrand;
    nr = (int*)malloc(sizeof(int)*n_rand);
    filenamesrand = (char***)malloc(sizeof(char**)*n_rand);
    for(i=0;i<n_rand;i++){
	nr[i] = Lines_ReadFile(Idxr_names[i],filenamesrand+i,NULL);
    }

    int nq;
    char **queries;
    nq = Lines_ReadFile(Idxq,&queries,NULL);

/*Read state sequences*/
    /*database*/
    char name_fullpath[FileNameLength_MAX];
    unsigned char **sseq;
    unsigned short **sseq_num;
    int *lf;
    int t;
    sseq = (unsigned char**)malloc(sizeof(unsigned char*)*nf);
    sseq_num = (unsigned short**)malloc(sizeof(unsigned short*)*nf);
    lf = (int*)malloc(sizeof(int)*nf);
    for(t=0;t<nf;t++){
	name_fullpath[0] = '\0';
	strcat(name_fullpath,Pf);
	strcat(name_fullpath,filenames[t]);
	lf[t] = Sseq_ReadFile(name_fullpath,sseq+t,sseq_num+t,NULL);
    }
    /*randomized baseline*/
    unsigned char ***sseq_rand;
    unsigned short ***sseq_num_rand;
    int **lr;
    sseq_rand = (unsigned char***)malloc(sizeof(unsigned char**)*n_rand);
    sseq_num_rand = (unsigned short***)malloc(sizeof(unsigned short**)*n_rand);
    lr = (int**)malloc(sizeof(int*)*n_rand);
    for(i=0;i<n_rand;i++){
	sseq_rand[i] = (unsigned char**)malloc(sizeof(unsigned char*)*nr[i]);
	sseq_num_rand[i] = (unsigned short**)malloc(sizeof(unsigned short*)*nr[i]);
	lr[i] = (int*)malloc(sizeof(int)*nr[i]);
	for(t=0;t<nr[i];t++){
	    name_fullpath[0] = '\0';
	    strcat(name_fullpath,Pr_names[i]);
	    strcat(name_fullpath,filenamesrand[i][t]);
	    lr[i][t] = Sseq_ReadFile(name_fullpath,sseq_rand[i]+t,sseq_num_rand[i]+t,NULL);
	}
    }
    /*query*/
    unsigned char **sseq_query;
    unsigned short **sseq_num_query;
    int *lq;
    sseq_query = (unsigned char**)malloc(sizeof(unsigned char*)*nq);
    sseq_num_query = (unsigned short**)malloc(sizeof(unsigned short*)*nq);
    lq = (int*)malloc(sizeof(int)*nq);
    for(t=0;t<nq;t++){
	name_fullpath[0] = '\0';
	strcat(name_fullpath,Pq);
	strcat(name_fullpath,queries[t]);
	lq[t] = Sseq_ReadFile(name_fullpath,sseq_query+t,sseq_num_query+t,NULL);
    }
/*Compute Cumsum of each chromesome*/
    int **sseq_cumnum;
    sseq_cumnum = (int**)malloc(sizeof(int*)*nf);
    for(t=0;t<nf;t++){
	sseq_cumnum[t] = (int*)malloc(sizeof(int)*(lf[t]+1));
	sseq_cumnum[t][0] = 0;
	for(i=0;i<lf[t];i++){
	    sseq_cumnum[t][i+1] = sseq_cumnum[t][i] + sseq_num[t][i];
	}
    }
/*Compute Cumsum of each query sequence*/
    int **sseq_query_cumnum;
    sseq_query_cumnum = (int**)malloc(sizeof(int*)*nq);
    for(t=0;t<nq;t++){
	sseq_query_cumnum[t] = (int*)malloc(sizeof(int)*(lq[t]+1));
	sseq_query_cumnum[t][0] = 0;
	for(i=0;i<lq[t];i++){
	    sseq_query_cumnum[t][i+1] = sseq_query_cumnum[t][i] + sseq_num_query[t][i];
	}
    }
/*Do the alignment for each query*/
    int h,hi;
    float **map_score;
    unsigned char **map_trace;
    float **map_score_rand;
    unsigned char **map_trace_rand;
    float **maxline;
    float **maxline_trace;
    float ***maxline_rand;
    maxline = (float**)malloc(sizeof(float*)*nf);
    maxline_trace = (float**)malloc(sizeof(float*)*nf);
    maxline_rand = (float***)malloc(sizeof(float**)*n_rand);
    for(i=0;i<n_rand;i++){
	maxline_rand[i] = (float**)malloc(sizeof(float*)*nr[i]);
    }
    unsigned char *sseq1;
    unsigned char *sseq2;
    unsigned short *sseq1_num;
    unsigned short *sseq2_num;
    int l1,l2;
    void *opt;
    float alpha;

    for(t=0;t<nq;t++){
	/*Open new file to record results*/
	char name_fullpath_region[FileNameLength_MAX];
	name_fullpath_region[0] = '\0';
	strcat(name_fullpath_region,Po);
	strcat(name_fullpath_region,queries[t]);
	strcat(name_fullpath_region,".algn");
	FILE *out_region = fopen(name_fullpath_region,"w");
	/*Align to random control*/
	float score_best_rand; /*score_best_rand is the average over all the randomized baselines*/
	score_best_rand = 0.0;
	for(hi=0;hi<n_rand;hi++){
	    for(h=0;h<nr[hi];h++){
		/*Deep copy of sseq1 and sseq2*/
		l1 = lq[t];
		sseq1 = (unsigned char*)malloc(sizeof(unsigned char)*l1);
		sseq1_num = (unsigned short*)malloc(sizeof(unsigned short)*l1);
		for(i=0;i<l1;i++){
		    sseq1[i] = sseq_query[t][i];
		    sseq1_num[i] = sseq_num_query[t][i];
		}
		l2 = lr[hi][h];
		sseq2 = (unsigned char*)malloc(sizeof(unsigned char)*l2);
		sseq2_num = (unsigned short*)malloc(sizeof(unsigned short)*l2);
		for(i=0;i<l2;i++){
		    sseq2[i] = sseq_rand[hi][h][i];
		    sseq2_num[i] = sseq_num_rand[hi][h][i];
		}
		/*Do alignment*/
		float *sseq1_num_align;
		float *sseq2_num_align;
		Custom_Init(&sseq1,&sseq1_num,&sseq1_num_align,&l1,&sseq2,&sseq2_num,&sseq2_num_align,&l2,&alpha,&opt,argv[3]);
		int l1p1,l2p1;
		l1p1 = l1+1;
		l2p1 = l2+1;
		map_score_rand = (float**)malloc(sizeof(float*)*l2p1);
		map_trace_rand = (unsigned char**)malloc(sizeof(unsigned char*)*l2p1);
		for(i=0;i<l2p1;i++){
		    map_score_rand[i] = (float*)malloc(sizeof(float)*l1p1);
		    map_trace_rand[i] = (unsigned char*)malloc(sizeof(unsigned char)*l1p1);
		}
		SWA_Even(sseq1, sseq2, l1, l2, Custom_MatchingFunction, Custom_GapFunction, alpha, NULL, opt, map_score_rand, map_trace_rand);
		/*Record the alignment score*/
		maxline_rand[hi][h] = (float*)malloc(sizeof(float)*l2p1);
		for(i=0;i<l2p1;i++){
		    maxline_rand[hi][h][i] = 0.0;
		    for(j=0;j<l1p1;j++){
			if(maxline_rand[hi][h][i]<map_score_rand[i][j]){
			    maxline_rand[hi][h][i] = map_score_rand[i][j];
			}
		    }
		}
		/*Free assigned space*/
		for(i=0;i<l2p1;i++){
		    free(map_score_rand[i]);
		    free(map_trace_rand[i]);
		}
		free(map_score_rand);
		free(map_trace_rand);
		Custom_Free(opt);
		free(sseq1);
		free(sseq1_num);
		free(sseq1_num_align);
		free(sseq2);
		free(sseq2_num);
		free(sseq2_num_align);
	    }
	
	    /*Compute best score align to randomized database*/
	    float score_best_rand_this;
	    score_best_rand_this = 0.0;
	    for(i=0;i<nr[hi];i++){
		for(j=0;j<lr[hi][i]+1;j++){
		    if(score_best_rand_this<maxline_rand[hi][i][j]){
			score_best_rand_this = maxline_rand[hi][i][j];
		    }
		}
	    }
	    for(h=0;h<nr[hi];h++){
		free(maxline_rand[hi][h]);
	    }
	    score_best_rand += score_best_rand_this;
	}
	score_best_rand = score_best_rand / n_rand;
	float th_this;
	th_this = l1*th_peak + score_best_rand;
	/*Align to database*/
	for(h=0;h<nf;h++){
	    /*Deep copy of sseq1 and sseq2*/
	    l1 = lq[t];
	    sseq1 = (unsigned char*)malloc(sizeof(unsigned char)*l1);
	    sseq1_num = (unsigned short*)malloc(sizeof(unsigned short)*l1);
	    for(i=0;i<l1;i++){
		sseq1[i] = sseq_query[t][i];
		sseq1_num[i] = sseq_num_query[t][i];
	    }
	    l2 = lf[h];
	    sseq2 = (unsigned char*)malloc(sizeof(unsigned char)*l2);
	    sseq2_num = (unsigned short*)malloc(sizeof(unsigned short)*l2);
	    for(i=0;i<l2;i++){
		sseq2[i] = sseq[h][i];
		sseq2_num[i] = sseq_num[h][i];
	    }
	    /*Do alignment*/
	    float *sseq1_num_align;
	    float *sseq2_num_align;
	    Custom_Init(&sseq1,&sseq1_num,&sseq1_num_align,&l1,&sseq2,&sseq2_num,&sseq2_num_align,&l2,&alpha,&opt,argv[3]);
	    int l1p1,l2p1;
	    l1p1 = l1+1;
	    l2p1 = l2+1;
	    map_score = (float**)malloc(sizeof(float*)*l2p1);
	    map_trace = (unsigned char**)malloc(sizeof(unsigned char*)*l2p1);
	    for(i=0;i<l2p1;i++){
		map_score[i] = (float*)malloc(sizeof(float)*l1p1);
		map_trace[i] = (unsigned char*)malloc(sizeof(unsigned char)*l1p1);
	    }	    
	    SWA_Even(sseq1, sseq2, l1, l2, Custom_MatchingFunction, Custom_GapFunction, alpha, NULL, opt, map_score, map_trace);
	    /*Record alignment score*/
	    maxline[h] = (float*)malloc(sizeof(float)*l2p1);
	    maxline_trace[h] = (float*)malloc(sizeof(float)*l2p1);
	    for(i=0;i<l2p1;i++){
		maxline[h][i] = 0.0;
		maxline_trace[h][i] = 0.0;
		for(j=0;j<l1p1;j++){
		    if(maxline[h][i]<map_score[i][j]){
			maxline[h][i] = map_score[i][j];
			maxline_trace[h][i] = map_score[i][j];
		    }
		}
	    }
	    /*Do traceback and Output high-score alignments*/
	    int pos_chr_this;
	    int pos_l1;
	    int chr_nothing;
	    struct pair_node align_this;
	    while(1){
		Find_Peak(maxline_trace+h,1,lf+h,&chr_nothing,&pos_chr_this);
		if(maxline_trace[h][pos_chr_this]>th_this){
		    /*print this alignment*/
		    pos_l1 = 0;
		    for(i=0;i<=l1;i++){
			if(map_score[pos_chr_this][i]>map_score[pos_chr_this][pos_l1]){
			    pos_l1 = i;
			}
		    }
		    Trace_Even(map_trace,pos_l1,pos_chr_this,&align_this);
		    fprintf(out_region,"chr:%d pos:%d-%d\nquery:%d-%d\nscore:%f\n",h,sseq_cumnum[h][align_this.next->p2],sseq_cumnum[h][pos_chr_this]-1,sseq_query_cumnum[t][align_this.next->p1],sseq_query_cumnum[t][pos_l1]-1,maxline_trace[h][pos_chr_this]-score_best_rand);
		    Print_Alignment_Sseq_Even(&align_this,sseq1,sseq[h],sseq1_num,sseq_num[h],out_region);
		    Free_Alignment(&align_this);
		}
		else{
		    break;
		}
		Remove_Peak(maxline_trace+h,lf+h,chr_nothing,pos_chr_this,l1,dist_peak_l,dist_peak_u);
	    }
	    /*Free assigned space*/
	    for(i=0;i<l2p1;i++){
		free(map_score[i]);
		free(map_trace[i]);
	    }
	    free(map_score);
	    free(map_trace);
	    Custom_Free(opt);
	    free(sseq1);
	    free(sseq1_num);
	    free(sseq1_num_align);
	    free(sseq2);
	    free(sseq2_num);
	    free(sseq2_num_align);
	}
	/*Output basic statistics of this query sequence*/
	fprintf(stdout,"%d\t%f",l1,score_best_rand);
	int chr;
	int pos;
	for(i=0;i<record_peak;i++){
	    Find_Peak(maxline,nf,lf,&chr,&pos);
	    fprintf(stdout,"\t%f",maxline[chr][pos]);
	    Remove_Peak(maxline,lf,chr,pos,l1,dist_peak_l,dist_peak_u);
	}
	fprintf(stdout,"\n");
	/*free maxline*/
	for(h=0;h<nf;h++){
	    free(maxline[h]);
	    free(maxline_trace[h]);
	}
	fclose(out_region);
    }
/*free assigned space*/
    for(i=0;i<n_rand;i++){
	free(Idxr_names[i]);
	free(Pr_names[i]);
    }
    free(Idxr_names);
    free(Pr_names);
    free(Idxr);
    free(Pr);
    return 1;
}




