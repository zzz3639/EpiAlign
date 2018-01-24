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
    u = SWF_MIN(p+l2*dist_peak_u,L[h]-1);
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

/*compute alignment score for two frequency vectors*/
float score_baseline(float *seqvec1, float *seqvec2, int StateMax)
{
    float ans;
    int i;
    ans = 0.0;
    for(i=0;i<StateMax;i++){
	ans += (seqvec1[i]-seqvec2[i])*(seqvec1[i]-seqvec2[i]);
    }
    ans = -sqrt(ans);
    return ans;
}

void malloc_seqvecarray(int w, int StateMax, float ***seqvecarray_p)
{
    seqvecarray_p[0] = (float**)malloc(sizeof(float*)*w);
    int i;
    for(i=0;i<w;i++){
	seqvecarray_p[0][i] = (float*)malloc(sizeof(float)*StateMax);
    }
    return;
}

void free_seqvecarray(int w, float **seqvecarray)
{
    int i;
    for(i=0;i<w;i++){
	free(seqvecarray[i]);
    }
    free(seqvecarray);
    return;
}

/*Compute frequency vector for a given sequence*/
void SeqFrequencyCount(unsigned char *seq, int m, float *seqvec, int StateMax)
{
    int i;
    for(i=0;i<StateMax;i++){
	seqvec[i] = 0.0;
    }
    float invm;
    invm = 1.0/m;
    for(i=0;i<m;i++){
	seqvec[seq[i]] += invm;
    }
    return;
}

/*compute seqvecarray for every segments of sseq, seqvecarray should have been malloced outside*/
/*if fit_model=="normal", then state frequency is counted on original sequence and compared, if fit_model=="compressed", then state frequency is counted after compression*/
void SeqVecPrepare(unsigned char *sseq, unsigned short *sseq_num, int l, int w, int k, float **seqvecarray, int StateMax, char *fit_model)
{
    unsigned char *seq;
    int m;
    m = Sseq2Seq(sseq,sseq_num,l,&seq);
    int ww;
    ww = w/k;
    int h;
    h = SWF_MAX(m/k+1-ww,0);
    int i,j;
    int s,t;
    unsigned char *seqtemp;
    seqtemp = (unsigned char*)malloc(sizeof(unsigned char)*w);
    unsigned char *sseqtemp;
    unsigned short *sseqtemp_num;
    int ltemp;
    
    for(i=0;i<h;i++){
	s = i*k;
	t = s+w;
	for(j=s;j<t;j++){
	    seqtemp[j-s] = seq[j];
	}
	if(strcmp(fit_model,"normal")==0){
            SeqFrequencyCount(seqtemp, w, seqvecarray[i], StateMax);
	}
	else if(strcmp(fit_model,"compressed")==0){
	    ltemp = Seq2Sseq(seqtemp,w,&sseqtemp,&sseqtemp_num);
	    SeqFrequencyCount(sseqtemp, ltemp, seqvecarray[i], StateMax);
	    free(sseqtemp);
	    free(sseqtemp_num);
	}
	else{
	    return;
	}
    }
    free(seqtemp);
    free(seq);
    return;
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
	printf("\nUsage:\n./run Path_Search_baseline Para_Search_baseline\n");
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
    int i,j,k;
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
/*Do fitting for each region*/
    int *sq;
    sq = (int*)malloc(sizeof(int)*nq);
    for(i=0;i<nq;i++){
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
    int **sr;
    sr = (int**)malloc(sizeof(int*)*n_rand);
    for(i=0;i<n_rand;i++){
	sr[i] = (int*)malloc(sizeof(int)*nr[i]);
    }
    float ***score_rand;
    score_rand = (float***)malloc(sizeof(float**)*n_rand);
    for(i=0;i<n_rand;i++){
	score_rand[i] = (float**)malloc(sizeof(float*)*nr[i]);
    }
    for(j=0;j<n_rand;j++){
	for(i=0;i<nr[j];i++){
	    sr[j][i] = SWF_MAX(cr[j][i]/wk+1-ww,0);
	    score_rand[j][i] = (float*)malloc(sizeof(float)*sr[j][i]);
	}
    }
    /*malloc frequency vector containers*/
    float ***seqvecarray_query;
    float ***seqvecarray_database;
    float ****seqvecarray_rand;
    seqvecarray_query = (float***)malloc(sizeof(float**)*nq);
    seqvecarray_database = (float***)malloc(sizeof(float**)*nf);
    seqvecarray_rand = (float****)malloc(sizeof(float***)*n_rand);
    for(i=0;i<n_rand;i++){
	seqvecarray_rand[i] = (float***)malloc(sizeof(float**)*nr[i]);
    }
    for(i=0;i<nq;i++){
	malloc_seqvecarray(sq[i], StateMax, seqvecarray_query+i);
    }
    for(i=0;i<nf;i++){
	malloc_seqvecarray(sf[i], StateMax, seqvecarray_database+i);
    }
    for(i=0;i<n_rand;i++){
	for(j=0;j<nr[i];j++){
	    malloc_seqvecarray(sr[i][j], StateMax, seqvecarray_rand[i]+j);
	}
    }
    /*prepare the frequency vectors*/
    for(i=0;i<nq;i++){
	SeqVecPrepare(sseq_query[i], sseq_num_query[i], lq[i], wm, wk, seqvecarray_query[i], StateMax, fit_model);
    }
    for(i=0;i<nf;i++){
	SeqVecPrepare(sseq[i], sseq_num[i], lf[i], wm, wk, seqvecarray_database[i], StateMax, fit_model);
    }
    for(i=0;i<n_rand;i++){
	for(j=0;j<nr[i];j++){
	    SeqVecPrepare(sseq_rand[i][j], sseq_num_rand[i][j], lr[i][j], wm, wk, seqvecarray_rand[i][j], StateMax, fit_model);
	}
    }
    /*iteration begins*/
    int tc;
    int r,rc;
    int in_rand;
    int aq,bq;
    int ar,br;
    int ad,bd;
    float maxscore_rand;
    float *maxscore_data;
    maxscore_data = (float*)malloc(sizeof(float)*record_peak);
    int chr_rand, pos_rand;
    int chr,pos;
    FILE *ScoreFile;
    FILE *AnnoFile;
    char annonum[FileNameLength_MAX];
    /*iterate through queries*/
    for(t=0;t<nq;t++){
	printf("\n\nRunning, this is query seq %d\n\n",t);
	name_fullpath[0] = '\0';
	strcat(name_fullpath,Po);
	strcat(name_fullpath,queries[t]);
	strcat(name_fullpath,".baselinescore");
	ScoreFile = fopen(name_fullpath,"w");
	for(tc=0;tc<sq[t];tc++){
	    aq = tc*wk;
	    bq = aq+wm;
	    /*align to randomized baseline*/
	    maxscore_rand = 0.0;
	    for(in_rand=0;in_rand<n_rand;in_rand++){
		for(r=0;r<nr[in_rand];r++){
		    for(rc=0;rc<sr[in_rand][r];rc++){
			score_rand[in_rand][r][rc] = score_baseline(seqvecarray_query[t][tc], seqvecarray_rand[in_rand][r][rc], StateMax);
		    }
		}
		Find_Peak(score_rand[in_rand], nr[in_rand], sr[in_rand], &chr_rand, &pos_rand);
		maxscore_rand += score_rand[in_rand][chr_rand][pos_rand];
	    }
	    maxscore_rand = maxscore_rand/n_rand;
	    /*align to database*/
	    for(r=0;r<nf;r++){
		for(rc=0;rc<sf[r];rc++){
		    score_data[r][rc] = score_baseline(seqvecarray_query[t][tc], seqvecarray_database[r][rc], StateMax);
		}
	    }
	    /*find peaks and output to files*/	
	    name_fullpath[0] = '\0';
	    strcat(name_fullpath,Po);
	    strcat(name_fullpath,queries[t]);
	    strcat(name_fullpath,"/");
	    strcat(name_fullpath,"a");
	    sprintf(annonum,"%d",aq);
	    strcat(name_fullpath,annonum);
	    strcat(name_fullpath,"b");
	    sprintf(annonum,"%d",bq);
	    strcat(name_fullpath,annonum);
	    strcat(name_fullpath,".sseq.anno");
	    AnnoFile = fopen(name_fullpath,"w");
            fprintf(AnnoFile,"Query region: ");
	    fprintf(AnnoFile,"chr%s, a%db%d\n\n",chrnames_query[t],aq,bq);

            for(i=0;i<record_peak;i++){
		Find_Peak(score_data, nf, sf, &chr, &pos);
		maxscore_data[i] = score_data[chr][pos];
		Remove_Peak(score_data, sf, chr, pos, ww, dist_peak_l, dist_peak_u);
		fprintf(AnnoFile,"Match%d: chr%s a%db%d, score: %f\n",i,chrnames_data[chr],pos*wk,pos*wk+wm,maxscore_data[i]);
	    }
            fprintf(ScoreFile,"%f",maxscore_rand);
	    for(i=0;i<record_peak;i++){
		fprintf(ScoreFile,"\t%f",maxscore_data[i]);
	    }
	    fprintf(ScoreFile,"\n");
            fclose(AnnoFile);

	}
	fclose(ScoreFile);
    }

/*free assigned space*/
    /*free filenames*/
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
    /*free sseq*/
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
    /*free cumnum*/
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
    /*free frequency vectors*/
    for(i=0;i<nq;i++){
	free_seqvecarray(sq[i],seqvecarray_query[i]);
    }
    free(seqvecarray_query);
    for(i=0;i<nf;i++){
	free_seqvecarray(sf[i],seqvecarray_database[i]);
    }
    free(seqvecarray_database);
    for(i=0;i<n_rand;i++){
	for(j=0;j<nr[i];j++){
	    free_seqvecarray(sr[i][j],seqvecarray_rand[i][j]);
	}
	free(seqvecarray_rand[i]);
    }
    free(seqvecarray_rand);
    /*free segment number vector*/
    free(sf);
    free(sq);
    for(i=0;i<n_rand;i++){
	free(sr[i]);
    }
    free(sr);
    /*free alignment score containers*/
    for(i=0;i<n_rand;i++){
	for(j=0;j<nr[i];j++){
	    free(score_rand[i][j]);
	}
    }
    for(i=0;i<n_rand;i++){
	free(score_rand[i]);
    }
    free(score_rand);
    /*free other spaces*/
    free(maxscore_data);
    free(Idxr);
    free(Pr);
    for(i=0;i<n_rand;i++){
	free(Idxr_names[i]);
	free(Pr_names[i]);
    }
    free(Idxr_names);
    free(Pr_names);
    free(nr);
    return 1;
}




