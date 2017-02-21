#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"EpiBLASTConstant.h"
#include"StateIO.h"

#define MatchKeyLen 3
#define ChrLetterMax 2

struct StatePair
{
    char c1;
    char c2;
    int l1;
    int l2;
};

struct algn
{
    float score;
    int chr;
    int s;
    int t;
    int qs;
    int qt;
};

void SwapFloat(float *x, float *y)
{
    float temp;
    temp = *x;
    *x = *y;
    *y = temp;
    return;
}

void SwapInt(int *x, int *y)
{
    int temp;
    temp = *x;
    *x = *y;
    *y = temp;
    return;
}


int ChoosePivot(int i,int j )
{
   return ((i+j)/2);
}

void QuickSort(float *A, int *idx, int m, int n)
{
    float key;
    int key_idx;
    int i,j,k;
    if(m<n){
	k = ChoosePivot(m,n);
	SwapFloat(A+m,A+k);
	SwapInt(idx+m,idx+k);
	key = A[m];
	i = m+1;
	j = n;
	while(i<=j){
	    while((i<=n)&&(A[i]<=key))
		i++;
	    while((j>=m)&&(A[j]>key))
		j--;
	    if(i<j){
		SwapFloat(A+i,A+j);
		SwapInt(idx+i,idx+j);
	    }
	}
	SwapFloat(A+m,A+j);
	SwapInt(idx+m,idx+j);
	QuickSort(A,idx,m,j-1);
	QuickSort(A,idx,j+1,n);
    }
    return;
}

void JumpEmpty(FILE *fp, char *emptychar)
{
    char c;
    c = fgetc(fp);
    while(emptychar[c]){
	c = fgetc(fp);
    }
    if(c!=EOF){
	fseek(fp,-1,SEEK_CUR);
    }
    return;
}

int JumpMatchKey(FILE *fp, char *matchkey)
{
    int i;
    int l=0;
    char c;
    char cp;
    cp='\n';
    c = fgetc(fp);
    i=0;
    while(c!=EOF){
	if(cp=='\n'&&c!='\n'){
	    l++;
	}
	cp = c;
	if(c==matchkey[i]){
	    i++;
	    if(i==MatchKeyLen){
		break;
	    }
	}
	else{
	    i=0;
	    if(c==matchkey[i]){
		i++;
	    }
	}
	c = fgetc(fp);
    }
    if(c!=EOF){
	fseek(fp,-MatchKeyLen,SEEK_CUR);
	l--;
    }
    return l;
}

int main(int argc, char **argv)
{
    if(argc<2){
	printf("\nUsage:\n");
	printf("./run AnnotationPara Prefix InputFile\n");
	return 1;
    }
    /*No two letters in matchkey should be identical in matchkey*/
    char matchkey[MatchKeyLen];
    matchkey[0]='c'; matchkey[1]='h'; matchkey[2]='r';
    int chrnamekeylen=3;
    char chrnamekey[chrnamekeylen+1];
    chrnamekey[0]='c'; chrnamekey[1]='h'; chrnamekey[2]='r'; chrnamekey[3]='\0';
    char emptychar[ASCII_MAX];
    int i,j;
    for(i=0;i<ASCII_MAX;i++){
	emptychar[i] = 0;
    }
    emptychar[' ']=1; emptychar['\n']=1; emptychar['\t']=1;

/*Read annotation parameters*/
    int qs;
    int qt;
    int lstate;
    int dbadd;
    int qadd;
    int chrquery;
    char ChrFile[FileNameLength_MAX];
    FILE *AnnotationFile;
    AnnotationFile = fopen(argv[1],"r");
    fscanf(AnnotationFile,"%d",&qs);
    fscanf(AnnotationFile,"%d",&qt);
    fscanf(AnnotationFile,"%d",&lstate);
    fscanf(AnnotationFile,"%d",&chrquery);
    fscanf(AnnotationFile,"%d",&qadd);
    fscanf(AnnotationFile,"%d",&dbadd);
    fscanf(AnnotationFile,"%s",ChrFile);
    fclose(AnnotationFile);
    char **chrnames;
    int nchr;
    nchr = Lines_ReadFile(ChrFile,&chrnames,NULL);
    char **chrletters;
    chrletters = (char**)malloc(sizeof(char*)*nchr);
    for(i=0;i<nchr;i++){
	chrletters[i] = (char*)malloc(sizeof(char)*(ChrLetterMax+1));
    }
    char *strsearchtemp;
    for(i=0;i<nchr;i++){
	strsearchtemp = strstr(chrnames[i],chrnamekey);
	j=chrnamekeylen;
	while(strsearchtemp[j]!='_'){
	    chrletters[i][j-chrnamekeylen] = strsearchtemp[j];
	    j++;
	}
	chrletters[i][j-chrnamekeylen] = '\0';
    }

/*Read input file*/
    int n=0;
    FILE *Input;
    char c;
    Input = fopen(argv[3],"r");
    c = fgetc(Input);
    i=0;
    while(c!=EOF){
	if(c==matchkey[i]){
	    i++;
	    if(i==MatchKeyLen){
		n++;
		i=0;
	    }
	}
	else{
	    i=0;
	    if(c==matchkey[i]){
		i++;
	    }
	}
	c = fgetc(Input);
    }
    fclose(Input);
    Input = fopen(argv[3],"r");
    struct algn *algnarray;
    algnarray = (struct algn*)malloc(sizeof(struct algn)*n);
    int k;
    float scorethis;
    int *lengthalgn;
    lengthalgn = (int*)malloc(sizeof(int)*n);
    for(i=0;i<n;i++){
	if(i>0){
	    lengthalgn[i-1] = JumpMatchKey(Input,matchkey);
	}
	fscanf(Input,"chr:%d",&k);
	algnarray[i].chr = k;
	JumpEmpty(Input,emptychar);
	fscanf(Input,"pos:%d-%d",&j,&k);
	algnarray[i].s = j;
	algnarray[i].t = k;
	JumpEmpty(Input,emptychar);
	fscanf(Input,"query:%d-%d",&j,&k);
	algnarray[i].qs = j;
	algnarray[i].qt = k;
	JumpEmpty(Input,emptychar);
	fscanf(Input,"score:%f",&scorethis);
	algnarray[i].score = scorethis;
    }
    lengthalgn[n-1] = JumpMatchKey(Input,matchkey);
    fclose(Input);
    /*Read Alignments*/
    struct StatePair **algncontent;
    algncontent = (struct StatePair**)malloc(sizeof(struct StatePair*)*n);
    for(i=0;i<n;i++){
	algncontent[i] = (struct StatePair*)malloc(sizeof(struct StatePair)*lengthalgn[i]);
    }
    Input = fopen(argv[3],"r");
    int t;
    for(i=0;i<n;i++){
	JumpMatchKey(Input,matchkey);
	fscanf(Input,"chr:%d",&k);
	JumpEmpty(Input,emptychar);
	fscanf(Input,"pos:%d-%d",&t,&k);
	JumpEmpty(Input,emptychar);
	fscanf(Input,"query:%d-%d",&t,&k);
	JumpEmpty(Input,emptychar);
	fscanf(Input,"score:%f",&scorethis);
	JumpEmpty(Input,emptychar);
        for(j=0;j<lengthalgn[i];j++){
	    fscanf(Input,"%c,%c\t%d,%d",&(algncontent[i][j].c1),&(algncontent[i][j].c2),&(algncontent[i][j].l1),&(algncontent[i][j].l2));
	    JumpEmpty(Input,emptychar);
	}
    }
    fclose(Input);
    /*printf("\n");
    for(i=0;i<n;i++){
	printf("%d,%d,%d,%d,%d,%f\n",algnarray[i].chr,algnarray[i].s,algnarray[i].t,algnarray[i].qs,algnarray[i].qt,algnarray[i].score);
    }*/
    /*printf("\n");
    for(i=0;i<n;i++){
	printf("%d\n",lengthalgn[i]);
    }*/
/*Convert to annotation file*/
    FILE *PrefixFile;
    char PrefixNewline;
    PrefixFile = fopen(argv[2],"r");
    c = fgetc(PrefixFile);
    while(c!=EOF){
	c = fgetc(PrefixFile);
    }
    fseek(PrefixFile,-1,SEEK_CUR);
    c = fgetc(PrefixFile);
    if(c=='\n'){
	PrefixNewline = 1;
    }
    else{
	PrefixNewline = 0;
    }
    fclose(PrefixFile);
    PrefixFile = fopen(argv[2],"r");
    while((c=fgetc(PrefixFile))!=EOF){
	fputc(c,stdout);
    }
    if(!PrefixNewline){
	fputc('\n',stdout);
    }
    fclose(PrefixFile);

    float *scorearray;
    int *idxarray;
    scorearray = (float*)malloc(sizeof(float)*n);
    idxarray = (int*)malloc(sizeof(int)*n);
    for(i=0;i<n;i++){
	scorearray[i] = algnarray[i].score;
	idxarray[i] = i;
    }
    QuickSort(scorearray,idxarray,0,n-1);

    fprintf(stdout,"Query: chr%s:%d-%d",chrletters[chrquery],qs*lstate+qadd,(qt+1)*lstate+qadd-1);
    fprintf(stdout,"\n%d matches found\n",n);
    for(i=0;i<n;i++){
	fprintf(stdout,"\nAlignment Pair:%d, score:%f\n",i,algnarray[idxarray[n-i-1]].score);
	fprintf(stdout,"Query:chr%s:%d-%d ",chrletters[chrquery],algnarray[idxarray[n-i-1]].qs*lstate+qadd,(algnarray[idxarray[n-i-1]].qt+1)*lstate+qadd-1);
	fprintf(stdout,"Database:chr%s:%d-%d\n",chrletters[algnarray[idxarray[n-i-1]].chr],algnarray[idxarray[n-i-1]].s*lstate+dbadd,(algnarray[idxarray[n-i-1]].t+1)*lstate+dbadd-1);
    }
    fprintf(stdout,"\nHere's the alignments:\n");
    int t1,t2;
    for(i=0;i<n;i++){
	fprintf(stdout,"\nAlignment Pair:%d, score:%f\n",i,algnarray[idxarray[n-i-1]].score);
	fprintf(stdout,"Query:chr%s:%d-%d ",chrletters[chrquery],algnarray[idxarray[n-i-1]].qs*lstate+qadd,(algnarray[idxarray[n-i-1]].qt+1)*lstate+qadd-1);
	fprintf(stdout,"Database:chr%s:%d-%d\n",chrletters[algnarray[idxarray[n-i-1]].chr],algnarray[idxarray[n-i-1]].s*lstate+dbadd,(algnarray[idxarray[n-i-1]].t+1)*lstate+dbadd-1);
	t1 = 0;
	t2 = 0;
	for(j=0;j<lengthalgn[idxarray[n-i-1]];j++){
	    fprintf(stdout,"%c,%c\t%d,%d",algncontent[idxarray[n-i-1]][j].c1,algncontent[idxarray[n-i-1]][j].c2,algncontent[idxarray[n-i-1]][j].l1*lstate,algncontent[idxarray[n-i-1]][j].l2*lstate);
	    fprintf(stdout,"\t%d-%d",(algnarray[idxarray[n-i-1]].qs+t1)*lstate+qadd, (algnarray[idxarray[n-i-1]].qs+t1+algncontent[idxarray[n-i-1]][j].l1)*lstate+qadd-1);
	    fprintf(stdout,",%d-%d\n",(algnarray[idxarray[n-i-1]].s+t2)*lstate+dbadd,(algnarray[idxarray[n-i-1]].s+t2+algncontent[idxarray[n-i-1]][j].l2)*lstate+dbadd-1);
	    t1 += algncontent[idxarray[n-i-1]][j].l1;
	    t2 += algncontent[idxarray[n-i-1]][j].l2;
	}
    }
    free(lengthalgn);
    for(i=0;i<nchr;i++){
	free(chrletters[i]);
    }
    free(chrletters);
    return 1;
}



