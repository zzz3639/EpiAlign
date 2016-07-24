#include"StateIO.h"

void MakeKey_Number(struct StateIO_Opt *opt)
{
    char c;
    const char cmax=ASCII_MAX-1;
    for(c='0';;c++){
	opt->key[c] = 0;
	if(c==cmax)
	    break;
    }
    for(c='0';c<='9';c++){
	opt->key[c] = 1;
    }
    return;
}

int FILE_CountLine(const char *FileName, struct StateIO_Opt *opt)
{
    FILE *fp=fopen(FileName,"r");
    int ans=0;
    char s;
    if(opt==NULL){
	s=0;
    }
    else{
	s=opt->s;
    }
    if(s==0){
	char c;
	char cp;
	c = fgetc(fp);
	cp = '\n';
	while(c!=EOF){
	    if(cp=='\n'){
		ans+=1;
	    }
	    cp = c;
	    c = fgetc(fp);
	}
    }
    else{
	char c;
	char cp;
	c = fgetc(fp);
	cp = '\n';
	while(c!=EOF){
	    if(cp=='\n'&&opt->key[c]){
		ans+=1;
	    }
	    cp = c;
	    c = fgetc(fp);
	}
    }
    fclose(fp);
    return ans;
}

int Sseq_ReadFile(const char *FileName, unsigned char** sseq, unsigned short** sseq_num, struct StateIO_Opt *opt)
{
    int n;
    struct StateIO_Opt Opt_Temp;
    if(opt==NULL){
	opt = &Opt_Temp;
	opt->s = 1;
	MakeKey_Number(opt);
	opt->n = FILE_CountLine(FileName,opt);
	opt->s = 0;
    }
    if(opt->n<1){
	struct StateIO_Opt opt2;
	opt2.s = 1;
	MakeKey_Number(&opt2);
	n = FILE_CountLine(FileName,&opt2);
    }
    else{
	n = opt->n;
    }
    unsigned char *sseq_temp;
    int *sseq_num_temp;
    sseq_temp = (unsigned char *)malloc(sizeof(unsigned char)*n);
    sseq_num_temp = (int *)malloc(sizeof(int)*n);
    FILE *fp=fopen(FileName,"r");
    char c;
    char cp;
    cp = '\n';
    c = fgetc(fp);
    int i=0;
    int k;
    while(c!=EOF){
	if(cp=='\n'&&opt->key[c]){
	    fseek(fp,-1,SEEK_CUR);
	    fscanf(fp,"%d",&k);
            sseq_temp[i] = k;
	    fscanf(fp,"%d",&k);
	    sseq_num_temp[i] = k;
	    i += 1;
	    c='\n';
	}
        if(i==n)
	    break;
	cp = c;
	c = fgetc(fp);
    }
    fclose(fp);
    int m=0;
    for(i=0;i<n;i++){
	k = sseq_num_temp[i];
	while(k>USHRT_MAX){
	    k -= USHRT_MAX;
	    m += 1;
	}
	m += 1;
    }
    if(opt->s==0){
	(*sseq) = (unsigned char *)malloc(sizeof(unsigned char)*m);
	(*sseq_num) = (unsigned short *)malloc(sizeof(unsigned short)*m);
    }
    k=0;
    for(i=0;i<n;i++){
	while(sseq_num_temp[i]>USHRT_MAX){
	    sseq_num_temp[i] -= USHRT_MAX;
	    *((*sseq)+k) = sseq_temp[i];
	    *((*sseq_num)+k) = USHRT_MAX;
	    k += 1;
	}
	*((*sseq)+k) = sseq_temp[i];
	*((*sseq_num)+k) = sseq_num_temp[i];
	k += 1;
    }
    free(sseq_temp);
    free(sseq_num_temp);
    return m;
}

int Seq_ReadFile(const char *FileName, unsigned char** seq, struct StateIO_Opt *opt)
{
    int n;
    struct StateIO_Opt Opt_Temp;
    if(opt==NULL){
	opt = &Opt_Temp;
	opt->s = 1;
	MakeKey_Number(opt);
	opt->n = FILE_CountLine(FileName,opt);
	opt->s = 0;
    }
    if(opt->n<1){
	struct StateIO_Opt opt2;
	opt2.s = 1;
	MakeKey_Number(&opt2);
	n = FILE_CountLine(FileName,&opt2);
    }
    else{
	n = opt->n;
    }
    if(opt->s==0){
	(*seq) = (unsigned char *)malloc(sizeof(unsigned char)*n);
    }
    FILE *fp=fopen(FileName,"r");
    char c;
    char cp;
    cp = '\n';
    c = fgetc(fp);
    int i=0;
    int k;
    while(c!=EOF){
	if(cp=='\n'&&opt->key[c]){
	    fseek(fp,-1,SEEK_CUR);
	    fscanf(fp,"%d",&k);
            *((*seq)+i) = k;
	    i += 1;
	    c='\n';
	}
	cp = c;
	c = fgetc(fp);
    }
    fclose(fp);
    return n;
}


int Lines_ReadFile(const char *InFileName, char ***Lines, struct StateIO_Opt *opt)
{
    int nr;
    if(opt==NULL){
	nr = FILE_CountLine(InFileName,NULL);
	(*Lines) = (char**)malloc(sizeof(char*)*nr);
    }
    else{
	/*ToBeFinished*/
    }
    int i=0;
    int j;
    FILE *fp = fopen(InFileName,"r");
    char c;
    char cp;
    cp = '\n';
    c = fgetc(fp);
    int m;
    while(c!=EOF){
	if(cp=='\n'&&c!='\n'){
	    /*count the length of this line*/
	    m = 0;
	    while(c!=EOF&&c!='\n'){
		m += 1;
                c = fgetc(fp);
	    }
	    if(c=='\n')
		fseek(fp,-m-1,SEEK_CUR);
	    else
		fseek(fp,-m,SEEK_CUR);
	    /*store this line*/
	    *((*Lines)+i) = (char*)malloc(sizeof(char)*(m+1));
	    for(j=0;j<m;j++){
		*(*((*Lines)+i)+j) = fgetc(fp);
	    }
	    *(*((*Lines)+i)+m) = '\0';
	    cp = *(*((*Lines)+i)+m-1);
	    c = fgetc(fp);
            i += 1;
	}
	else{
	    cp = c;
	    c = fgetc(fp);
	}
    }
    fclose(fp);
    return nr;
}




