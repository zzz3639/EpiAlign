#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<limits.h>
#include<string.h>
#include"WatermanFun.h"
#include"StateIO.h"

int main(int argc, char **argv)
{
    if(argc<1){
	printf("\n./run parameters\n");
        printf("Parameters:\nm\nk\nlist path\ndatabase path\noutput folder path\n");
        return 1;
    }
/*read parameters*/
    int m;
    int k;
    char Pf[FileNameLength_MAX];
    char Pl[FileNameLength_MAX];
    char Po[FileNameLength_MAX];
    FILE *para=fopen(argv[1],"r");
    fscanf(para,"%d",&m);
    fscanf(para,"%d",&k);
    fscanf(para,"%s",Pl);
    fscanf(para,"%s",Pf);
    fscanf(para,"%s",Po);
    fclose(para);
/*read filenames*/
    int nf;
    char **filenames;
    nf = Lines_ReadFile(Pl,&filenames,NULL);

/*Read state sequences*/
    int i,j;
    int t;
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

/*Compute Cumsum of each chromesome*/
    unsigned short **sseq_cumnum;
    sseq_cumnum = (unsigned short**)malloc(sizeof(unsigned short*)*nf);
    for(t=0;t<nf;t++){
	sseq_cumnum[t] = (unsigned short*)malloc(sizeof(unsigned short)*(lf[t]+1));
	sseq_cumnum[t][0] = 0;
	for(i=0;i<lf[t];i++){
	    sseq_cumnum[t][i+1] = sseq_cumnum[t][i] + sseq_num[t][i];
	}
    }
/*Do cut on each short sequences defined by m and k*/
    unsigned char *seq;
    unsigned char *seq2;
    seq2 = (unsigned char*)malloc(sizeof(unsigned char)*m);
    int h;
    unsigned char *sseq2;
    unsigned short *sseq2_num;
    int l2;
    for(t=0;t<nf;t++){
	seq = (unsigned char*)malloc(sizeof(unsigned char)*lfull[t]);
	h=0;
	for(i=0;i<lf[t];i++){
	    for(j=h;j<h+sseq_num[t][i];j++){
		seq[j] = sseq[t][i];
	    }
	    h = h+sseq_num[t][i];
	}
	int a,b;
        FILE *SseqIndex;
        char index_fullpath[FileNameLength_MAX];
	index_fullpath[0] = '\0';
	strcat(index_fullpath,Po);
	strcat(index_fullpath,filenames[t]);
	strcat(index_fullpath,".idx");
        SseqIndex = fopen(index_fullpath,"w");
	a = 0;
	b = m;
	while(b<lfull[t]){
	    for(i=0;i<m;i++){
		seq2[i] = seq[a+i];
	    }
	    l2 = Seq2Sseq(seq2,b-a,&sseq2,&sseq2_num);
	    

	    char name_fullpath_region[FileNameLength_MAX];
	    char string_temp[FileNameLength_MAX];
	    char name_suffix[FileNameLength_MAX];
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

            /*Print file name to index*/
	    strcpy(name_suffix,"a");
            sprintf(string_temp,"%d",a);
	    strcat(name_suffix,string_temp);
	    strcat(name_suffix,"b");
	    sprintf(string_temp,"%d",b);
	    strcat(name_suffix,string_temp);
            strcat(name_suffix,".sseq");
            fprintf(SseqIndex,"%s\n",name_suffix);
            /*Output to sseq file*/
	    strcpy(name_suffix,name_fullpath_region);
            strcat(name_suffix,".sseq");
	    FILE *out_region = fopen(name_suffix,"w");
	    for(i=0;i<l2;i++){
		fprintf(out_region,"%d %d\n",(sseq2[i]),sseq2_num[i]);
	    }
	    fclose(out_region);
            /*Output to seq file*/
	    strcpy(name_suffix,name_fullpath_region);
            strcat(name_suffix,".seq");
	    out_region = fopen(name_suffix,"w");
	    for(i=0;i<m;i++){
		fprintf(out_region,"%d\n",(seq2[i]));
	    }
	    fclose(out_region);

            a += k;
            b += k;
	}
        fclose(SseqIndex);
	free(seq);
    }

    
/*free spaces*/
    for(i=0;i<nf;i++){
	free(filenames[i]);
    }
    free(filenames);
    free(lf);
    free(lfull);
    for(i=0;i<nf;i++){
	free(sseq[i]);
        free(sseq_num[i]);
        free(sseq_cumnum[i]);
    }
    free(sseq);
    free(sseq_num);
    free(sseq_cumnum);
    return 1;
}






