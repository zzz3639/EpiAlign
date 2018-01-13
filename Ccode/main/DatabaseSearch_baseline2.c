#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"WatermanFun.h"
#include"StateIO.h"
#include"EpiBLASTConstant.h"

#define ChrLetterMax 2
#define chrnamekeylen 3
/*in case of very long line for randomized baseline*/
#define n_rand_max 5

/*As scores are all negative, regions removed are assigned by -STATE_MAX*/
void Remove_Peak(float **maxline, int *L, int h, int p, int l2, int dist_peak_l, int dist_peak_u)
{
    int l;
    int u;
    l = SWF_MAX(p-l2*dist_peak_l,0);
    u = SWF_MIN(p+l2*dist_peak_u,L[h]);
    int i;
    for(i=l;i<=u;i++){
	maxline[h][i] = -1.0*STATE_MAX;
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
	for(j=0;j<L[i];j++){
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

int Num_Parsing(char *strin)
{
    int i,j;
    char c;
    i=0;
    j=1;
    while((c=strin[i])!='\0'){
	if(c==' '){
	    j++;
	}
	i++;
    }
    return j;
}

void Ending_Parsing(char *strin)
{
    int i;
    i=0;
    while(strin[i]!='\0'){
        if(strin[i]=='\n'){
	    strin[i] = '\0';
	    break;
	}
	i++;
    }
    return;
}

void Path_Parsing(char *strin, char **strarrayout)
{
    char c;
    int i,j,k;
    i=0; j=0; k=0;
    while((c=strin[i])!='\0'){
	strarrayout[j][i-k] = c;
	if(c==' '){
	    strarrayout[j][i-k] = '\0';
	    k = i+1;
	    j++;
	}
	i++;
    }
    strarrayout[j][i-k] = '\0';
    return;
}

int main(int argc, char **argv)
{
    if(argc<3){
	printf("\nUsage:\n./run Path_Search_baseline Para_Search_baseline Para_Annotation_baseline\n");
	return 1;
    }

/*read searching parameters from file*/
    FILE *para_file=fopen(argv[2],"r");
    int wm;
    int wk;
    int dist_peak_l;
    int dist_peak_u;
    int record_peak;
    char fit_model[FileNameLength_MAX];
    fscanf(para_file,"%d",&wm);
    fscanf(para_file,"%d",&wk);
    fscanf(para_file,"%d",&dist_peak_l);
    fscanf(para_file,"%d",&dist_peak_u);
    fscanf(para_file,"%d",&record_peak);
    fscanf(para_file,"%s",fit_model);
    fclose(para_file);
    int ww=wm/wk;
/*read annotation parameters*/
    /*
    FILE *para_anno=fopen(argv[3],"r");
    fclose(para_anno);*/
/*read paths*/
    char Idxq[FileNameLength_MAX];
    char Pq[FileNameLength_MAX];
    char Idxf[FileNameLength_MAX];
    char Pf[FileNameLength_MAX];
    char *Idxr;
    char *Pr;
    char Po[FileNameLength_MAX];
    FILE *seq_paths=fopen(argv[1],"r");
    Idxr = (char*)malloc(sizeof(char)*FileNameLength_MAX*n_rand_max);
    Pr = (char*)malloc(sizeof(char)*FileNameLength_MAX*n_rand_max);
    fscanf(seq_paths,"%s",Idxq);
    fscanf(seq_paths,"%s",Pq);
    fscanf(seq_paths,"%s",Idxf);
    fscanf(seq_paths,"%s",Pf);
    fgetc(seq_paths);
    fgets(Idxr,FileNameLength_MAX*n_rand_max,seq_paths);
    fgets(Pr,FileNameLength_MAX*n_rand_max,seq_paths);
    fscanf(seq_paths,"%s",Po);
    fclose(seq_paths);

/*Parsing the randomized database paths and indexes*/
    int n_rand;
    int i,j;
    n_rand = Num_Parsing(Idxr);
    Ending_Parsing(Idxr);
    Ending_Parsing(Pr);
    char **Idxr_names;
    char **Pr_names;
    Idxr_names = (char**)malloc(sizeof(char*)*n_rand);
    Pr_names = (char**)malloc(sizeof(char*)*n_rand);
    for(i=0;i<n_rand;i++){
	Idxr_names[i] = (char*)malloc(sizeof(char)*FileNameLength_MAX);
	Pr_names[i] = (char*)malloc(sizeof(char)*FileNameLength_MAX);
    }
    Path_Parsing(Idxr,Idxr_names);
    Path_Parsing(Pr,Pr_names);

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

/*find chromosome letters for query and database*/
    char chrnamekey[chrnamekeylen+1];
    chrnamekey[0]='c'; chrnamekey[1]='h'; chrnamekey[2]='r'; chrnamekey[3]='\0';
    char **chrnames_query;
    char **chrnames_data;
    chrnames_query = (char**)malloc(sizeof(char*)*nq);
    for(i=0;i<nq;i++){
	chrnames_query[i] = (char*)malloc(sizeof(char)*(ChrLetterMax+1));
    }
    chrnames_data = (char**)malloc(sizeof(char*)*nf);
    for(i=0;i<nf;i++){
	chrnames_data[i] = (char*)malloc(sizeof(char)*(ChrLetterMax+1));
    }
    char *strsearchtemp;
    for(i=0;i<nq;i++){
	strsearchtemp = strstr(queries[i],chrnamekey);
	j=chrnamekeylen;
	while(strsearchtemp[j]!='_'){
	    chrnames_query[i][j-chrnamekeylen] = strsearchtemp[j];
	    j++;
	}
	chrnames_query[i][j-chrnamekeylen] = '\0';
    }
    for(i=0;i<nf;i++){
	strsearchtemp = strstr(filenames[i],chrnamekey);
	j=chrnamekeylen;
	while(strsearchtemp[j]!='_'){
	    chrnames_data[i][j-chrnamekeylen] = strsearchtemp[j];
	    j++;
	}
	chrnames_data[i][j-chrnamekeylen] = '\0';
    }

/*Count StateMax*/
    int StateMax;
    char StateMaxArray[STATE_MAX];
    for(i=0;i<STATE_MAX;i++){
	StateMaxArray[i] = 0;
    }
    for(t=0;t<nq;t++){
	for(i=0;i<lq[t];i++){
	    StateMaxArray[sseq_query[t][i]]=1;
	}
    }
    for(t=0;t<nf;t++){
	for(i=0;i<lf[t];i++){
	    StateMaxArray[sseq[t][i]]=1;
	}
    }
    for(j=0;j<n_rand;j++){
	for(t=0;t<nr[j];t++){
	    for(i=0;i<lr[j][t];i++){
		StateMaxArray[sseq_rand[j][t][i]]=1;
	    }
	}
    }
    StateMax = 0;
    for(i=0;i<STATE_MAX;i++){
	if(StateMaxArray[i]==1){
	    StateMax = i+1;
	}
    }
/*Compute length of each chromesome*/
    int *cf;
    cf = (int*)malloc(sizeof(int)*nf);
    for(t=0;t<nf;t++){
	cf[t] = 0;
	for(i=0;i<lf[t];i++){
	    cf[t] += sseq_num[t][i];
	}
    }
    int *cq;
    cq = (int*)malloc(sizeof(int)*nq);
    for(t=0;t<nq;t++){
	cq[t] = 0;
	for(i=0;i<lq[t];i++){
	    cq[t] += sseq_num_query[t][i];
	}
    }
    int **cr;
    cr = (int**)malloc(sizeof(int*)*n_rand);
    for(i=0;i<n_rand;i++){
	cr[i] = (int*)malloc(sizeof(int)*nr[i]);
    }
    for(j=0;j<n_rand;j++){
	for(t=0;t<nr[j];t++){
	    cr[j][t] = 0;
	    for(i=0;i<lr[j][t];i++){
		cr[j][t] += sseq_num_rand[j][t][i];
	    }
	}
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
/*Do fitting for each region*/ /*we are here*/
    int *sq;
    sq = (int*)malloc(sizeof(int)*nq);
    float **score_query;
    score_query = (float**)malloc(sizeof(float*)*nq);
    for(i=0;i<nq;i++){
	score_query[i] = (float*)malloc(sizeof(float)*SWF_MAX(cq[i]/wk+1-ww,0));
	sq[i] = SWF_MAX(cq[i]/wk+1-ww,0);
    }
    int *sf;
    sf = (int*)malloc(sizeof(int)*nf);
    float **score_data;
    score_data = (float**)malloc(sizeof(float*)*nf);
    for(i=0;i<nf;i++){
	score_data[i] = (float*)malloc(sizeof(float)*SWF_MAX(cf[i]/wk+1-ww,0));
	sf[i] = SWF_MAX(cf[i]/wk+1-ww,0);
    }
    int *sr;
    sr = (int*)malloc(sizeof(int)*nr);
    float **score_rand;
    score_rand = (float**)malloc(sizeof(float*)*nr);
    for(i=0;i<nr;i++){
	score_rand[i] = (float*)malloc(sizeof(float)*SWF_MAX(cr[i]/wk+1-ww,0));
	sr[i] = SWF_MAX(cr[i]/wk+1-ww,0);
    }
    int k;
    int r;
    int aoq,boq;
    int aiq,biq;
    int aor,bor;
    int air,bir;
    int aod,bod;
    int aid,bid;
    unsigned char *seq_q_temp;
    unsigned char *seq_r_temp;
    unsigned char *seq_db_temp;
    unsigned char *seq_q_this;
    unsigned char *seq_r_this;
    unsigned char *seq_db_this;
    seq_q_this = (unsigned char*)malloc(sizeof(unsigned char)*wm);
    seq_r_this = (unsigned char*)malloc(sizeof(unsigned char)*wm);
    seq_db_this = (unsigned char*)malloc(sizeof(unsigned char)*wm);
    float maxscore_rand;
    float *maxscore_data;
    maxscore_data = (float*)malloc(sizeof(float)*record_peak);
    FILE *ScoreFile;
    FILE *AnnoFile;
    int chr;
    int pos;
    int chr_rand;
    int pos_rand;
    unsigned char *sseq_q_this;
    unsigned short *sseq_num_q_this;
    int l_q_this;
    unsigned char *sseq_db_this;
    unsigned short *sseq_num_db_this;
    int l_db_this;
    char annonum[FileNameLength_MAX];
    for(t=0;t<nq;t++){
        printf("\n\nRunning, this is query seq %d\n\n",t);
	name_fullpath[0] = '\0';
	strcat(name_fullpath,Po);
	strcat(name_fullpath,queries[t]);
	strcat(name_fullpath,".baselinescore");
        ScoreFile = fopen(name_fullpath,"w");
        aoq = 0;
        boq = wm;
        aiq = aoq/wk;
        biq = boq/wk;
	seq_q_temp = (unsigned char*)malloc(sizeof(unsigned char)*cq[t]);
	k=0;
	for(i=0;i<lq[t];i++){
	    for(j=0;j<sseq_num_query[t][i];j++){
		seq_q_temp[k] = sseq_query[t][i];
		k++;
	    }
	}
        while(boq<=cq[t]){
            printf("Segment %d\n",aiq);
	    /*initialize*/
	    /*find state sequence of this region*/
            for(i=aoq;i<boq;i++){
                seq_q_this[i-aoq] = seq_q_temp[i];
	    }
	    /*do fitting for this region on randomized dataset*/
	    for(r=0;r<nr;r++){
		seq_r_temp = (unsigned char*)malloc(sizeof(unsigned char)*cr[r]);
                k=0;
		for(i=0;i<lr[r];i++){
		    for(j=0;j<sseq_num_rand[r][i];j++){
			seq_r_temp[k] = sseq_rand[r][i];
			k++;
		    }
		}
		aor = 0;
		bor = wm;
		air = aor/wk;
		bir = bor/wk;
		while(bor<=cr[r]){
		    for(i=aor;i<bor;i++){
			seq_r_this[i-aor] = seq_r_temp[i];
		    }
		    score_rand[r][air] = score_baseline(seq_r_this,wm,seq_q_this,wm,fit_model,StateMax);
		    /*printf("\n%d\t%d\t%d\n",r,cr[r],air);*/
		    aor += wk;
		    bor += wk;
		    air += 1;
		    bir += 1;
		}
		free(seq_r_temp);
	    }
	    /*find best score for randomized sequence*/
            Find_Peak(score_rand, nr, sr, &chr_rand, &pos_rand);
            maxscore_rand = score_rand[chr_rand][pos_rand];
	    /*do fitting for this region on true dataset*/
            for(r=0;r<nf;r++){
		seq_db_temp = (unsigned char*)malloc(sizeof(unsigned char)*cf[r]);
		k=0;
		for(i=0;i<lf[r];i++){
                    for(j=0;j<sseq_num[r][i];j++){
			seq_db_temp[k] = sseq[r][i];
			k++;
		    }
		}
		aod = 0;
		bod = wm;
		aid = aod/wk;
		bid = bod/wk;
		while(bod<=cf[r]){
		    for(i=aod;i<bod;i++){
			seq_db_this[i-aod] = seq_db_temp[i];
		    }
		    score_data[r][aid] = score_baseline(seq_db_this,wm,seq_q_this,wm,fit_model,StateMax);
                    aod += wk;
		    bod += wk;
		    aid += 1;
		    bid += 1;
		}
		free(seq_db_temp);
	    }
	    /*find peaks*/
	    name_fullpath[0] = '\0';
	    strcat(name_fullpath,Po);
	    strcat(name_fullpath,queries[t]);
	    strcat(name_fullpath,"/");
	    strcat(name_fullpath,"a");
	    sprintf(annonum,"%d",aoq);
	    strcat(name_fullpath,annonum);
	    strcat(name_fullpath,"b");
	    sprintf(annonum,"%d",boq);
	    strcat(name_fullpath,annonum);
	    strcat(name_fullpath,".sseq.anno");
	    AnnoFile = fopen(name_fullpath,"w");
            fprintf(AnnoFile,"Query region: ");
	    fprintf(AnnoFile,"chr%s, a%db%d\n\n",chrnames_query[t],aoq,boq);

            for(i=0;i<record_peak;i++){
		Find_Peak(score_data, nf, sf, &chr, &pos);
		maxscore_data[i] = score_data[chr][pos];
		Remove_Peak(score_data, sf, chr, pos, ww, dist_peak_l, dist_peak_u);
		fprintf(AnnoFile,"Match%d: chr%s a%db%d, score: %f\n",i,chrnames_data[chr],pos*wk,pos*wk+wm,maxscore_data[i]);
	    }
	    fclose(AnnoFile);
	    /*record results*/
	    fprintf(ScoreFile,"%f",maxscore_rand);
	    for(i=0;i<record_peak;i++){
		fprintf(ScoreFile,"\t%f",maxscore_data[i]);
	    }
	    fprintf(ScoreFile,"\n");
	    aoq += wk;
	    boq += wk;
	    aiq += 1;
	    biq += 1;
        }
	free(seq_q_temp);
	fclose(ScoreFile);
    }

/*free assigned space*/
    for(i=0;i<nf;i++){
	free(filenames[i]);
    }
    free(filenames);
    for(j=0;j<n_rand;j++){
	for(i=0;i<nr[j];i++){
	    free(filenamesrand[j][i]);
	}
    }
    for(i=0;i<n_rand;i++){
	free(filenamesrand[i]);
    }
    free(filenamesrand);
    for(i=0;i<nq;i++){
	free(queries[i]);
    }
    free(queries);
    for(i=0;i<nf;i++){
	free(sseq[i]);
	free(sseq_num[i]);
    }
    free(sseq);
    free(sseq_num);
    for(j=0;j<n_rand;j++){
	for(i=0;i<nr[j];i++){
	    free(sseq_rand[j][i]);
	    free(sseq_num_rand[j][i]);
	}
    }
    for(i=0;i<n_rand;i++){
	free(sseq_rand[i]);
	free(sseq_num_rand[i]);
    }
    free(sseq_rand);
    free(sseq_num_rand);
    for(i=0;i<nq;i++){
	free(sseq_query[i]);
	free(sseq_num_query[i]);
    }
    free(sseq_query);
    free(sseq_num_query);
    for(i=0;i<nf;i++){
	free(sseq_cumnum[i]);
    }
    free(sseq_cumnum);
    for(i=0;i<nq;i++){
	free(sseq_query_cumnum[i]);
    }
    free(sseq_query_cumnum);
    free(cf);
    free(cq);
    for(i=0;i<n_rand;i++){
	free(cr[i]);
    }
    free(cr);
    free(sf);
    free(sq);
    free(sr);
    free(Idxr);
    free(Pr);
    free(seq_q_this);
    free(seq_r_this);
    free(seq_db_this);
    free(maxscore_data);
    for(i=0;i<n_rand;i++){
	free(Idxr_names[i]);
	free(Pr_names[i]);
    }
    free(Idxr_names);
    free(Pr_names);
    return 1;
}




