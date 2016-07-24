#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<limits.h>
#include<string.h>
#include"WatermanFun.h"
#include"StateIO.h"

#define dist_peak_l 1
#define dist_peak_u 2

#define record_peak 3
#define th_peak 0.25

float MatchScore_Naive_This(unsigned char a, unsigned char b, int i, int j, void *opt)
{
    return (a==b)?2.0:-2.0;
}

void Remove_Peak(float **maxline, int n, int *L, int h, int p, int l2)
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
    if(argc<3){
	printf("\n./run filelist filelistrand parameters remapping\n");
        return 1;
    }
/*read parameters*/
    int m;
    float alpha;
    int k;
    char Pf[FileNameLength_MAX];
    char Pr[FileNameLength_MAX];
    char Po[FileNameLength_MAX];
    FILE *para=fopen(argv[3],"r");
    fscanf(para,"%d",&m);
    fscanf(para,"%d",&k);
    fscanf(para,"%f",&alpha);
    fscanf(para,"%s",Pf);
    fscanf(para,"%s",Pr);
    fscanf(para,"%s",Po);
    fclose(para);
/*read re-mapping dictionary*/
    int i,j;
    int t;
    int Dn;
    Dn = FILE_CountLine(argv[4],NULL);
    FILE *DicFile=fopen(argv[4],"r");
    unsigned char *DicState;
    DicState = (unsigned char*)malloc(sizeof(unsigned char)*(Dn+1));
    for(i=0;i<Dn;i++){
	fscanf(DicFile,"%d",&j);
	fscanf(DicFile,"%d",&t);
	DicState[j] = t;
    }
    fclose(DicFile);
/*read filenames*/
    int nf;
    char **filenames;
    nf = Lines_ReadFile(argv[1],&filenames,NULL);
    int nr;
    char **filenamesrand;
    nr = Lines_ReadFile(argv[2],&filenamesrand,NULL);

/*Read state sequences*/
    char name_fullpath[FileNameLength_MAX];
    unsigned char **sseq;
    unsigned short **sseq_num;
    int *lf;
    sseq = (unsigned char**)malloc(sizeof(unsigned char*)*nf);
    sseq_num = (unsigned short**)malloc(sizeof(unsigned short*)*nf);
    lf = (int*)malloc(sizeof(int)*nf);
    for(t=0;t<nf;t++){
	name_fullpath[0] = '\0';
	strcat(name_fullpath,Pf);
	strcat(name_fullpath,filenames[t]);
	lf[t] = Sseq_ReadFile(name_fullpath,sseq+t,sseq_num+t,NULL);
    }
    unsigned char **sseq_rand;
    unsigned short **sseq_num_rand;
    int *lr;
    sseq_rand = (unsigned char**)malloc(sizeof(unsigned char*)*nr);
    sseq_num_rand = (unsigned short**)malloc(sizeof(unsigned short*)*nr);
    lr = (int*)malloc(sizeof(int)*nr);
    for(t=0;t<nr;t++){
	name_fullpath[0] = '\0';
	strcat(name_fullpath,Pr);
	strcat(name_fullpath,filenamesrand[t]);
	lr[t] = Sseq_ReadFile(name_fullpath,sseq_rand+t,sseq_num_rand+t,NULL);
    }

/*print lenghts of the files*/
    int *lfull;
    lfull = (int*)malloc(sizeof(int)*nf);
    for(t=0;t<nf;t++){
	lfull[t] = 0;
	for(i=0;i<lf[t];i++){
	    lfull[t] += sseq_num[t][i];
	}
	printf("\n%d %d\n",t,lfull[t]);
    }
    int *lfull_rand;
    lfull_rand = (int*)malloc(sizeof(int)*nr);
    for(t=0;t<nr;t++){
	lfull_rand[t] = 0;
	for(i=0;i<lr[t];i++){
	    lfull_rand[t] += sseq_num_rand[t][i];
	}
    }

/*re-map this sequence*/
    unsigned char *seq_full_temp;
    for(t=0;t<nf;t++){
	for(i=0;i<lf[t];i++){
	    sseq[t][i] = DicState[sseq[t][i]];
	}
	Sseq2Seq(sseq[t],sseq_num[t],lf[t],&seq_full_temp);
	free(sseq[t]);
	free(sseq_num[t]);
	lf[t] = Seq2Sseq(seq_full_temp,lfull[t],sseq+t,sseq_num+t);
	free(seq_full_temp);
    }
    for(t=0;t<nr;t++){
	for(i=0;i<lr[t];i++){
	    sseq_rand[t][i] = DicState[sseq_rand[t][i]];
	}
	Sseq2Seq(sseq_rand[t],sseq_num_rand[t],lr[t],&seq_full_temp);
	free(sseq_rand[t]);
	free(sseq_num_rand[t]);
	lr[t] = Seq2Sseq(seq_full_temp,lfull_rand[t],sseq_rand+t,sseq_num_rand+t);
	free(seq_full_temp);
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
/*Do the alignment on each short sequences defined by m and k*/
    unsigned char *seq;
    unsigned char *seq2;
    seq2 = (unsigned char*)malloc(sizeof(unsigned char)*m);
    int h;
    float **map_score;
    unsigned char **map_trace;
    float **map_score_rand;
    unsigned char **map_trace_rand;
    unsigned char *sseq2;
    unsigned short *sseq2_num;
    int l2;
    float **maxline;
    float **maxline_trace;
    float **maxline_rand;
    maxline = (float**)malloc(sizeof(float*)*nf);
    maxline_trace = (float**)malloc(sizeof(float*)*nf);
    maxline_rand = (float**)malloc(sizeof(float*)*nr);
    for(i=0;i<nf;i++){
	maxline[i] = (float*)malloc(sizeof(float)*(lf[i]+1));
	maxline_trace[i] = (float*)malloc(sizeof(float)*(lf[i]+1));
    }
    for(i=0;i<nr;i++){
	maxline_rand[i] = (float*)malloc(sizeof(float)*(lr[i]+1));
    }
    FILE *out_p;
    for(t=0;t<nf;t++){
	name_fullpath[0] = '\0';
	strcat(name_fullpath,Po);
	strcat(name_fullpath,filenames[t]);
	strcat(name_fullpath,".score");
	out_p = fopen(name_fullpath,"w");
	seq = (unsigned char*)malloc(sizeof(unsigned char)*lfull[t]);
	h=0;
	for(i=0;i<lf[t];i++){
	    for(j=h;j<h+sseq_num[t][i];j++){
		seq[j] = sseq[t][i];
	    }
	    h = h+sseq_num[t][i];
	}
	int a,b;
	a = 0;
	b = m;
	while(b<lfull[t]){
	    for(i=0;i<m;i++){
		seq2[i] = seq[a+i];
	    }
	    l2 = Seq2Sseq(seq2,b-a,&sseq2,&sseq2_num);
	    
	    for(h=0;h<nr;h++){
		/*malloc the space*/
		map_score_rand = (float**)malloc(sizeof(float*)*(lr[h]+1));
		map_trace_rand = (unsigned char**)malloc(sizeof(unsigned char*)*(lr[h]+1));
		for(i=0;i<lr[h]+1;i++){
		    map_score_rand[i] = (float*)malloc(sizeof(float)*(l2+1));
		    map_trace_rand[i] = (unsigned char*)malloc(sizeof(unsigned char)*(l2+1));
		}
		/*Do alignment*/
		SWA_Even(sseq2, sseq_rand[h], l2, lr[h], MatchScore_Naive_This, GapScore_Naive, alpha, NULL, NULL, map_score_rand, map_trace_rand);
		/*compute maximum matching score*/
		for(i=0;i<lr[h]+1;i++){
		    maxline_rand[h][i] = 0.0;
		    for(j=0;j<l2+1;j++){
                        if(maxline_rand[h][i]<map_score_rand[i][j]){
			    maxline_rand[h][i] = map_score_rand[i][j];
			}
		    }
		}
		/*free the space*/
		for(i=0;i<lr[h]+1;i++){
		    free(map_score_rand[i]);
		    free(map_trace_rand[i]);
		}
		free(map_score_rand);
		free(map_trace_rand);
	    }
	    float score_best_rand;
	    score_best_rand = 0.0;
	    for(i=0;i<nr;i++){
		for(j=0;j<lr[i]+1;j++){
		    if(score_best_rand<maxline_rand[i][j]){
			score_best_rand = maxline_rand[i][j];
		    }
		}
	    }

	    char name_fullpath_region[FileNameLength_MAX];
	    char string_temp[FileNameLength_MAX];
	    name_fullpath_region[0] = '\0';
	    strcat(name_fullpath_region,Po);
	    strcat(name_fullpath_region,filenames[t]);
	    strcat(name_fullpath_region,"/");
	    strcat(name_fullpath_region,"a");
	    sprintf(string_temp,"%d",a);
	    strcat(name_fullpath_region,string_temp);
	    strcat(name_fullpath_region,"b");
	    sprintf(string_temp,"%d",b);
	    strcat(name_fullpath_region,string_temp);
	    FILE *out_region = fopen(name_fullpath_region,"w");
	    float th_this;
	    th_this = l2*th_peak + score_best_rand;
	    for(h=0;h<nf;h++){
                /*Malloc the space*/
		map_score = (float**)malloc(sizeof(float*)*(lf[h]+1));
		map_trace = (unsigned char**)malloc(sizeof(unsigned char*)*(lf[h]+1));
		for(i=0;i<lf[h]+1;i++){
		    map_score[i] = (float*)malloc(sizeof(float)*(l2+1));
		    map_trace[i] = (unsigned char*)malloc(sizeof(unsigned char)*(l2+1));
		}
		/*Do alignment*/
		SWA_Even(sseq2, sseq[h], l2, lf[h], MatchScore_Naive_This, GapScore_Naive, alpha, NULL, NULL, map_score, map_trace);
		/*compute maximum matching score*/
		for(i=0;i<lf[h]+1;i++){
		    maxline[h][i] = 0.0;
		    maxline_trace[h][i] = maxline[h][i];
		    for(j=0;j<l2+1;j++){
                        if(maxline[h][i]<map_score[i][j]){
			    maxline[h][i] = map_score[i][j];
			    maxline_trace[h][i] = maxline[h][i];
			}
		    }
		}
		/*Do trace back*/
		int pos_chr_this;
		int pos_l2;
		int chr_nothing;
		struct pair_node align_this;
		while(1){
		    Find_Peak(maxline_trace+h,1,lf+h,&chr_nothing,&pos_chr_this);
		    if(maxline_trace[h][pos_chr_this]>th_this){
			/*print this alignment*/
			pos_l2 = 0;
			for(i=0;i<=l2;i++){
			    if(map_score[pos_chr_this][i]>map_score[pos_chr_this][pos_l2]){
				pos_l2 = i;
			    }
			}
			fprintf(out_region,"chr:%d pos:%d\n score:%f\n",h,sseq_cumnum[h][pos_chr_this],maxline_trace[h][pos_chr_this]-score_best_rand);
			Trace_Even(map_trace,pos_l2,pos_chr_this,&align_this);
			Print_Alignment_Sseq_Even(&align_this,sseq2,sseq[h],sseq2_num,sseq_num[h],out_region);
			Free_Alignment(&align_this);
		    }
		    else{
			break;
		    }
		    Remove_Peak(maxline_trace+h,1,lf+h,chr_nothing,pos_chr_this,l2);
		}
		/*Free the space*/
		for(i=0;i<lf[h]+1;i++){
		    free(map_score[i]);
		    free(map_trace[i]);
		}
		free(map_score);
		free(map_trace);
	    }
	    fclose(out_region);
	    fprintf(out_p,"%d\t%f",l2,score_best_rand);
	    int chr;
	    int pos;
	    for(i=0;i<record_peak;i++){
		Find_Peak(maxline,nf,lf,&chr,&pos);
		fprintf(out_p,"\t%f",maxline[chr][pos]);
		Remove_Peak(maxline,nf,lf,chr,pos,l2);
	    }
	    fprintf(out_p,"\n");
            a += k;
            b += k;
	}
	free(seq);
	fclose(out_p);
    }

    
/*free spaces*/
    for(i=0;i<nf;i++){
	free(filenames[i]);
    }
    free(filenames);
    for(i=0;i<nf;i++){
	free(filenamesrand[i]);
    }
    free(filenamesrand);
    free(lr);
    free(lf);
    free(lfull);
    free(lfull_rand);
    for(i=0;i<nf;i++){
	free(maxline[i]);
	free(maxline_trace[i]);
    }
    free(maxline);
    free(maxline_trace);
    for(i=0;i<nr;i++){
	free(maxline_rand[i]);
    }
    free(maxline_rand);
    free(DicState);
    for(i=0;i<nf;i++){
	free(sseq[i]);
        free(sseq_num[i]);
        free(sseq_cumnum[i]);
    }
    free(sseq);
    free(sseq_num);
    free(sseq_cumnum);
    for(i=0;i<nr;i++){
	free(sseq_rand[i]);
	free(sseq_num_rand[i]);
    }
    free(sseq_rand);
    free(sseq_num_rand);
    return 1;
}






