#include"WatermanFun.h"

char State2Char(unsigned char c)
{
    return (c<27)?(char)(c)+'`':(char)(c)+'&';
}

unsigned char Char2State(char c)
{
    return (c>'`')?(unsigned char)c-(unsigned char)96:(unsigned char)c-(unsigned char)38;
}

float MatchScore_Naive(unsigned char a, unsigned char b, int i, int j, void* opt)
{
    return (a==b)?2.0:-2.0;
}

float GapScore_Naive(float alpha, int i, void *opt)
{
    return alpha;
}

/*Smith-Waterman function for even gap penalty*/
/*By default, seq1 is query and seq2 is database*/
void SWA_Even(unsigned char* seq1, unsigned char* seq2, int l1, int l2, MatchingFunction MF, GapFunction GF, float alpha, struct Map_State_Even* map_init, void* opt, float** map_score, unsigned char** map_trace)
{
    int i,j;
    int lp1=l1+1;
    int lp2=l2+1;
    /*Initialize*/
    for(i=0;i<lp1;i++)
	map_trace[0][i]=0;
    for(i=0;i<lp2;i++)
	map_trace[i][0]=0;
    if(map_init==NULL){
	for(i=0;i<lp1;i++)
	    map_score[0][i]=0;
	for(i=0;i<lp2;i++)
	    map_score[i][0]=0;
    }
    else{
	for(i=0;i<lp1;i++)
	    map_score[0][i]=map_init->score[i];
	for(i=0;i<lp2;i++)
	    map_score[i][0]=0;
    }
    /*Update the matrix*/
    for(i=1;i<lp2;i++){
        int iloop,jloop;
	iloop=i-1;
	for(j=1;j<lp1;j++){
	    jloop=j-1;
	    float matchthis=MF(seq1[jloop],seq2[iloop],jloop,iloop,opt);
	    float scoretemp;
	    map_trace[i][j]=0;
	    map_score[i][j]=0;
            /*match*/
	    scoretemp = map_score[iloop][jloop] + matchthis;
	    if(scoretemp>0.0){
		map_score[i][j] = scoretemp;
		map_trace[i][j] = 1;
	    }
	    /*deletion in sequence 1*/
	    scoretemp = map_score[i][jloop] - GF(alpha,jloop,opt);
	    if(scoretemp>map_score[i][j]){
		map_score[i][j] = scoretemp;
		map_trace[i][j] = 2;
	    }
	    /*deletion in sequence 2*/
	    scoretemp = map_score[iloop][j] - GF(alpha,iloop,opt);
	    if(scoretemp>map_score[i][j]){
		map_score[i][j] = scoretemp;
		map_trace[i][j] = 3;
	    }
	}
    }
    return;
}

int Seq2Sseq(unsigned char *seq, int n, unsigned char **sseq, unsigned short **sseq_len)
{
    if(n<1){
	sseq=NULL;
	sseq_len=NULL;
	return 0;
    }

    int i,j,k,m;
    unsigned char s;
    int t,r;
    int USHRT_m1=USHRT_MAX-1;
    /*Compute short sequence length*/
    k=0;
    s=seq[0];
    t=1;
    for(i=1;i<n;i++){
	if(seq[i]==s){
	    t+=1;
	    continue;
	}
	k+=(t+USHRT_m1)/USHRT_MAX;
	s=seq[i];
	t=1;
    }
    k+=(t+USHRT_m1)/USHRT_MAX;
    /*assign short sequence*/
    *(sseq)=(unsigned char *)malloc(sizeof(unsigned char)*k);
    *(sseq_len)=(unsigned short *)malloc(sizeof(unsigned short)*k);
    s=seq[0];
    t=1;
    j=0;
    for(i=1;i<n;i++){
	if(seq[i]==s){
	    t+=1;
	    continue;
	}
        /*assign short sequence here*/
	if(t<=USHRT_MAX){
	    *(*sseq+j)=s;
	    *(*sseq_len+j)=t;
	    j+=1;
	}
	else{
	    r=(t+USHRT_m1)/USHRT_MAX;
            for(m=0;m<r-1;m++){
		*(*sseq+j)=s;
		*(*sseq_len+j)=USHRT_MAX;
		j+=1;
	    }
	    *(*sseq+j)=s;
	    *(*sseq_len+j)=(t-1)%USHRT_MAX+1;
	    j+=1;
	}
        s=seq[i];
	t=1;
    }
    /*assign the last region*/
    if(t<=USHRT_MAX){
	*(*sseq+j)=s;
	*(*sseq_len+j)=t;
	j+=1;
    }
    else{
	r=(t+USHRT_m1)/USHRT_MAX;
	for(m=0;m<r-1;m++){
	    *(*sseq+j)=s;
	    *(*sseq_len+j)=USHRT_MAX;
	    j+=1;
	}
	*(*sseq+j)=s;
	*(*sseq_len+j)=(t-1)%USHRT_MAX+1;
	j+=1;
    }
    return k;
}

int Sseq2Seq(unsigned char *sseq, unsigned short *sseq_num, int n, unsigned char **seq)
{
    int l=0;
    int i,j;
    for(i=0;i<n;i++){
	l += sseq_num[i];
    }
    (*seq) = (unsigned char*)malloc(sizeof(unsigned char)*l);
    l=0;
    for(i=0;i<n;i++){
	for(j=0;j<sseq_num[i];j++){
	    *((*seq)+l) = sseq[i];
	    l++;
	}
    }
    return l;
}

void Malloc_Map_Compact(unsigned char* seq1, unsigned char* seq2, unsigned short* len1, unsigned short* len2, int l1, int l2, float*** map_score, unsigned char*** map_trace)
{
    int i,j,k;
    int lp1=l1+1;
    int lp2=l2+1;
    int n1=0,n2=0;
    for(i=0;i<l1;i++)
	n1+=len1[i];
    for(i=0;i<l2;i++)
	n2+=len2[i];
    *map_score=(float**)malloc(sizeof(float*)*(n2+1));
    *map_trace=(unsigned char**)malloc(sizeof(unsigned char*)*(n2+1));

    /*malloc the columns*/
    int np1=n1+1;
    char *M=(char*)malloc(sizeof(char)*(n2+1));
    for(i=0;i<=n2;i++)
	M[i]=0;
    **map_score=(float*)malloc(sizeof(float)*np1);
    **map_trace=(unsigned char*)malloc(sizeof(unsigned char)*np1);
    M[0]=1;
    k=0;
    for(i=0;i<l2;i++){
	k+=len2[i];
	*(*map_score+k)=(float*)malloc(sizeof(float)*np1);
	*(*map_trace+k)=(unsigned char *)malloc(sizeof(unsigned char)*np1);
        M[k]=1;
    }
    for(i=0;i<=n2;i++){
	if(M[i]==1)
	    continue;
	*(*map_score+i)=(float*)malloc(sizeof(float)*lp1);
	*(*map_trace+i)=(unsigned char*)malloc(sizeof(unsigned char)*lp1);
    }
    free(M);
    return;
}

void Free_Map(int n, float** map_score, unsigned char** map_trace)
{
    int i;
    for(i=0;i<n+1;i++){
	free(map_score[i]);
	free(map_trace[i]);
    }
    free(map_score);
    free(map_trace);
    return;
}
/*
void SWA_Compact_Even(unsigned char* seq1, unsigned char* seq2, unsigned short* len1, unsigned short* len2, int l1, int l2, MatchingFunction MF, float alpha, void *opt, float** map_score, unsigned char** map_trace)
{
    int i,j,k;
    int lp1=l1+1;
    int lp2=l2+1;
    int n1=0,n2=0;
    for(i=0;i<l1;i++)
	n1+=len1[i];
    for(i=0;i<l2;i++)
	n2+=len2[i];
    int np1=n1+1;
    int np2=n2+1;
    char *M;
    M=(char*)malloc(sizeof(char)*np2);
    for(i=0;i<np2;i++)
	M[i]=0;
    M[0]=1;
    k=0;
    for(i=0;i<l2;i++){
	k+=len2[i];
        M[k]=1;
    }
    int *Index1=(int*)malloc(sizeof(int)*(l1+1));
    Index1[0]=0;
    for(i=0;i<l1;i++)
	Index1[i+1]=Index1[i]+len1[i];
    for(i=0;i<np2;i++){
	if(M[k]==1){
	    for(j=0;j<np1;j++){
	        map_score[i][j]=0.0;
	        map_trace[i][j]=0;
	    }
	}
	else{
	    for(j=0;j<lp1;j++){
	        map_score[i][j]=0.0;
	        map_trace[i][j]=0;
	    }
	}
    }
    int im,jm;
    float score_m,score_d1,score_d2;
    unsigned short r1,r2;
    int t1,t2;
    r2=len2[0];
    t2=0;
    for(i=0;i<n2;i++){
	if(M[i]==0){
	    int iplusr2=i+r2;
	    for(j=0;j<l1;j++){
		r1=len1[j];
		score_d1=map_score[i][j] - r1*alpha;
		int jplus1=j+1;
		if(score_d1>map_score[i][jplus1]){
		    map_score[i][jplus1]=score_d1;
		    map_trace[i][jplus1]=2;
		}
		score_d2=map_score[i][j] - r2*alpha;
		if(score_d2>map_score[iplusr2][Index1[j]]){
		    map_score[iplusr2][Index1[j]]=score_d2;
		    map_trace[iplusr2][Index1[j]]=3;
		}
		if(r1>=r2){
		    score_m=map_score[i][j] + r2*MF(seq2[t2],seq1[j],opt);
		    int jplusr2=Index1[j]+r2;
                    if(score_m>map_score[iplusr2][jplusr2]){
			map_score[iplusr2][jplusr2]=score_m;
			map_trace[iplusr2][jplusr2]=1;
		    }
		}
		else{
		    score_m=map_score[i][j] + r1*MF(seq2[t2],seq1[j],opt);
		    int iplusr1=i+r1;
		    if(score_m>map_score[iplusr1][jplus1]){
			map_score[iplusr1][jplus1]=score_m;
			map_trace[iplusr1][jplus1]=1;
		    }
		}
	    }
	    score_d2=map_score[i][l1] - r2*alpha;
	    if(score_d2>map_score[iplusr2][n1]){
		map_score[iplusr2][n1]=score_d2;
		map_trace[iplusr2][n1]=3;
	    }
	}
	else{
	    t1=0;
	    r1=len1[t1];
	    int iplusr2=i+r2;
	    for(j=0;j<n1;j++){
		int jplusr1=j+r1;
		int jplus1=j+1;
		score_d1=map_score[i][j] - r1*alpha;
		if(score_d1>map_score[i][jplusr1]){
		    map_score[i][jplusr1]=score_d1;
		    map_trace[i][jplusr1]=2;
		}
		score_d2=map_score[i][j] - r2*alpha;
		if(score_d2>map_score[iplusr2][j]){
		    map_score[iplusr2][j]=score_d2;
		    map_trace[iplusr2][j]=3;
		}
		if(r1>=r2){
		    score_m=map_score[i][j] + r2*MF(seq2[t2],seq1[t1],opt);
		    int jplusr2=j+r2;
                    if(score_m>map_score[iplusr2][jplusr2]){
			map_score[iplusr2][jplusr2]=score_m;
			map_trace[iplusr2][jplusr2]=1;
		    }
		}
		else{
		    score_m=map_score[i][j] + r2*MF(seq2[t2],seq1[t1],opt);
		    int iplusr1=i+r1;
		    if(score_m>map_score[iplusr1][jplus1]){
			map_score[iplusr1][jplus1]=score_m;
			map_trace[iplusr1][jplus1]=1;
		    }
		}
		r1-=1;
		if(r1==0){
		    t1+=1;
		    r1=len1[t1];
		}
	    }
	    score_d2=map_score[i][n1] - r2*alpha;
	    if(score_d2>map_score[iplusr2][n1]){
		map_score[iplusr2][n1]=score_d2;
		map_trace[iplusr2][n1]=3;
	    }
	}
	r2-=1;
	if(r2==0){
	    t2+=1;
	    r2=len2[t2];
	}
    }
    free(M);
    free(Index1);
    return;
}
*/

void SWA_Compact_Even(unsigned char* seq1, unsigned char* seq2, unsigned short* len1, unsigned short* len2, int l1, int l2, MatchingFunction MF, float alpha, void *opt, float** map_score, unsigned char** map_trace)
{
    /*Define parameters*/
    int i,j,k;
    int lp1=l1+1;
    int lp2=l2+1;
    int n1=0,n2=0;
    for(i=0;i<l1;i++)
	n1+=len1[i];
    for(i=0;i<l2;i++)
	n2+=len2[i];
    int np1=n1+1;
    int np2=n2+1;
    /*Initialize*/
    char *M;
    M=(char*)malloc(sizeof(char)*np2);
    for(i=0;i<np2;i++)
	M[i]=0;
    M[0]=1;
    k=0;
    for(i=0;i<l2;i++){
	k+=len2[i];
        M[k]=1;
    }
    int *Index1=(int*)malloc(sizeof(int)*(l1+1));
    Index1[0]=0;
    for(i=0;i<l1;i++)
	Index1[i+1]=Index1[i]+len1[i];
    for(i=0;i<np2;i++){
	if(M[k]==1){
	    for(j=0;j<np1;j++){
	        map_score[i][j]=0.0;
	        map_trace[i][j]=0;
	    }
	}
	else{
	    for(j=0;j<lp1;j++){
	        map_score[i][j]=0.0;
	        map_trace[i][j]=0;
	    }
	}
    }
    /*Update the matrix*/
    int im,jm;
    float score_m,score_d1,score_d2;
    unsigned short r1,r2;
    int t1,t2;
    t2=0;
    r2=1;
    for(i=1;i<np2;i++){
	int iminus1=i-1;
	int iminusr2=i-r2;
	if(M[i]){
	    t1=0;
            r1=1;
	    for(j=1;j<np1;j++){
		int jminus1=j-1;
		int jminusr1=j-r1;
		int jminusr2=j-r2;
		int iminusr1=i-r1;
		score_d1 = map_score[i][jminus1] - alpha;
		if(r1>=r2){
		    score_m = map_score[iminusr2][jminusr2] + r2*MF(seq1[t1],seq2[t2],t1,t2,opt);
		}
		else{
		    score_m = map_score[iminusr1][t1] + r1*MF(seq1[t1],seq2[t2],t1,t2,opt);
		}
		if(r1==len1[t1]){
		    t1=t1+1;
		    score_d2 = map_score[iminus1][t1] - alpha;
		    r1=1;
		}
		else{
		    score_d2 = map_score[iminusr2][j] - r2*alpha;
		    r1+=1;
		}
		if(score_m>map_score[i][j]){
		    map_score[i][j]=score_m;
		    map_trace[i][j]=1;
		}
		if(score_d1>map_score[i][j]){
		    map_score[i][j]=score_d1;
		    map_trace[i][j]=2;
		}
		if(score_d2>map_score[i][j]){
		    map_score[i][j]=score_d2;
		    map_trace[i][j]=3;
		}
	    }
	}
	else{
	    for(j=1;j<lp1;j++){
		int jminus1=j-1;
		int jminusr2=Index1[j]-r2;
		int iminusr1=i-r1;
		r1=len1[jminus1];
		score_d1 = map_score[i][jminus1] - r1*alpha;
		if(r1>=r2){
		    score_m = map_score[iminusr2][jminusr2] + r2*MF(seq1[jminus1],seq2[t2],jminus1,t2,opt);
		}
		else{
		    score_m = map_score[iminusr1][jminus1] + r1*MF(seq1[jminus1],seq2[t2],jminus1,t2,opt);
		}
		if(M[iminus1]){
		    score_d2 = map_score[iminus1][Index1[j]] - alpha;
		}
		else{
		    score_d2 = map_score[iminus1][j] - alpha;
		}
		if(score_m>map_score[i][j]){
		    map_score[i][j]=score_m;
		    map_trace[i][j]=1;
		}
		if(score_d1>map_score[i][j]){
		    map_score[i][j]=score_d1;
		    map_trace[i][j]=2;
		}
		if(score_d2>map_score[i][j]){
		    map_score[i][j]=score_d2;
		    map_trace[i][j]=3;
		}
	    }
	}
	r2+=1;
	if(r2>len2[t2]){
	    t2+=1;
	    r2=1;
	}
    }
    free(M);
    free(Index1);
    return;
}


/*bow_seq is malloced inside this function*/
/*state 0 is reserved state here*/
/*This function does the same thing with function Seq2Bow, but time cost of this function is higher*/
int Seq2Bow(unsigned char *seq, int n, int m, int k, struct word_node** bow_seq)
{
    int i,j;
    int w=m/k;
    int a=n/k;
    int l=a-w+1;
    n = a*k;
    if(l<1){
	*bow_seq = NULL;
	return 0;
    }
    (*bow_seq) = (struct word_node*)malloc(sizeof(struct word_node)*l);
    /*Construct the length-k dictionary*/
    struct word_node_extend *dic;
    struct word_node_extend *dic_node;
    struct word_node_extend **dic_last;
    dic_last = (struct word_node_extend**)malloc(sizeof(struct word_node_extend*)*STATE_MAX);
    dic = (struct word_node_extend*)malloc(sizeof(struct word_node_extend)*STATE_MAX);
    for(i=0;i<STATE_MAX;i++){
	dic[i].next=NULL;
	dic[i].p=-1;
        dic_last[i]=dic+i;
    }
    unsigned char c;
    for(i=0;i<n;i++){
	c=seq[i];
	if(i - dic_last[c]->p * k < k){
	    dic_last[c]->n += 1.0;
	    continue;
	}
	dic_node = (struct word_node_extend*)malloc(sizeof(struct word_node_extend));
        dic_last[c]->next = dic_node;
	dic_node->next = NULL;
	dic_node->n = 1.0;
	dic_node->p = i/k;
	dic_last[c] = dic_node;
    }

    /*Compute the length-m bag of words*/
    struct word_node *word_node;
    struct word_node **word_last;
    word_last = (struct word_node**)malloc(sizeof(struct word_node*)*l);
    int p;
    int s,t;
    for(i=0;i<l;i++){
        word_last[i] = (*bow_seq) + i;
	word_last[i]->c = 0;
	word_last[i]->next = NULL;
    }
    for(i=1;i<STATE_MAX;i++){
	c = i;
	dic_node = dic[i].next;
	while(dic_node!=NULL){
	    p = dic_node->p;
	    s = SWF_MAX(p-w+1,0);
	    t = SWF_MIN(p+1,l);
	    for(j=s;j<t;j++){
		if(word_last[j]->c!=c){
		    word_node = (struct word_node *)malloc(sizeof(struct word_node));
		    word_node->c = c;
		    word_node->n = dic_node->n;
		    word_node->next = NULL;
		    word_last[j]->next = word_node;
		    word_last[j] = word_node;
		}
		else{
		    word_last[j]->n += dic_node->n;
		}
	    }
	    dic_node = dic_node->next;
	}
    }

    /*free the memory in this function*/
    for(i=0;i<STATE_MAX;i++){
	dic_last[i] = (dic+i)->next;
	while(dic_last[i]!=NULL){
	    dic_node = dic_last[i];
	    dic_last[i] = dic_node->next;
	    free(dic_node);
	}
    }
    free(dic);
    free(dic_last);
    free(word_last);
    return l;
}


/*bow_seq is malloced inside this function*/
/*state 0 is a reserved state in this function*/
int Seq2Bow_Slow(unsigned char *seq, int n, int m, int k, struct word_node** bow_seq)
{
    int i,j;
    int w=m/k;
    int t=n/k;
    int l=t-w+1;
    n = t*k;
    if(l<1){
	*bow_seq = NULL;
	return 0;
    }
    (*bow_seq) = (struct word_node*)malloc(sizeof(struct word_node)*l);
    /*Construct the dictionary*/
    struct word_node_extend *dic;
    struct word_node_extend *dic_node;
    struct word_node_extend *dic_last[STATE_MAX];
    dic = (struct word_node_extend*)malloc(sizeof(struct word_node_extend)*STATE_MAX);
    for(i=0;i<STATE_MAX;i++){
	dic[i].next=NULL;
	dic[i].p=-1;
        dic_last[i]=dic+i;
    }
    unsigned char c;
    for(i=0;i<n;i++){
	c=seq[i];
	if(i - dic_last[c]->p * k < k){
	    dic_last[c]->n += 1.0;
	    continue;
	}
	dic_node = (struct word_node_extend*)malloc(sizeof(struct word_node_extend));
        dic_last[c]->next = dic_node;
	dic_node->next = NULL;
	dic_node->n = 1.0;
	dic_node->p = i/k;
	dic_last[c] = dic_node;
    }
    /*compute length-k bag og word*/
    struct word_node *bow_seq_k;
    bow_seq_k = (struct word_node*)malloc(sizeof(struct word_node)*t);
    struct word_node **bow_last;
    bow_last = (struct word_node**)malloc(sizeof(struct word_node*)*t);
    struct word_node *word_pointer;
    for(i=0;i<t;i++){
	bow_seq_k[i].c = 0;
        bow_seq_k[i].next = NULL;
	bow_last[i] = bow_seq_k+i;
    }
    for(c=1;;c++){
	dic_node = dic[c].next;
	while(dic_node!=NULL){
	    i = dic_node->p;
	    word_pointer = (struct word_node*)malloc(sizeof(struct word_node));
	    word_pointer->c = c;
	    word_pointer->n = dic_node->n;
	    word_pointer->next = NULL;
	    bow_last[i]->next = word_pointer;
	    bow_last[i] = word_pointer;
	    dic_node = dic_node->next;
	}
        if(c==STATE_MAX-1){
	    break;
        }
    }

    /*compute length-m bag of word*/
    struct word_node *bow_seq_temp;
    bow_seq_temp = (struct word_node*)malloc(sizeof(struct word_node)*l);
    for(i=0;i<l;i++){
	(*bow_seq+i)->c = 0;
        (*bow_seq+i)->n = 0;
	(*bow_seq+i)->next = NULL;
	bow_seq_temp[i].next = NULL;
    }
    for(i=0;i<w;i++){ 
	Bow_Plus(*bow_seq,bow_seq_k+i,*bow_seq);
    }
    for(i=1;i<l;i++){
	Bow_Plus((*bow_seq+i-1),bow_seq_k+i+w-1,bow_seq_temp+i);
	Bow_Minus(bow_seq_temp+i,bow_seq_k+i-1,(*bow_seq+i));
    }

    /*free the memory in this function*/
    for(i=0;i<STATE_MAX;i++){
	dic_last[i] = (dic+i)->next;
	while(dic_last[i]!=NULL){
	    dic_node = dic_last[i];
	    dic_last[i] = dic_node->next;
	    free(dic_node);
	}
    }
    free(dic);
    for(i=0;i<t;i++){
	Free_Bow(bow_seq_k+i);
    }
    free(bow_seq_k);
    free(bow_last);
    for(i=0;i<l;i++){
	Free_Bow(bow_seq_temp+i);
    }
    free(bow_seq_temp);
    return l;
}


float MatchScore_Bow_Resemblance(struct word_node *bow1, struct word_node *bow2, void* opt)
{
    struct word_node *p1;
    struct word_node *p2;
    p1 = bow1->next;
    p2 = bow2->next;
    /*a is number of intersection and b is union*/
    float a,b;
    a=0.0;
    b=0.0;
    while(p1!=NULL&&p2!=NULL){
        if(p1->c==p2->c){
	    a += 1.0;
	    b += 1.0;
	    p1 = p1->next;
	    p2 = p2->next;
	    continue;
	}
	if(p1->c<p2->c){
	    b += 1.0;
	    p1 = p1->next;
	    continue;
	}
	if(p1->c>p2->c){
	    b += 1.0;
	    p2 = p2->next;
	    continue;
	}
    }
    if(p1!=NULL){
	while(p1!=NULL){
	    b += 1.0;
	    p1 = p1->next;
	}
    }
    if(p2!=NULL){
	while(p2!=NULL){
	    b += 1.0;
	    p2 = p2->next;
	}
    }
    float ans;
    ans=3.0*a/b-1.0;
    return ans;
}

float MatchScore_Bow_Resemblance_ByNumber(struct word_node *bow1, struct word_node *bow2, void* opt)
{
    struct word_node *p1;
    struct word_node *p2;
    p1 = bow1->next;
    p2 = bow2->next;
    /*a is number of intersection and b is union*/
    float a,b;
    a=0.0;
    b=0.0;
    while(p1!=NULL&&p2!=NULL){
	if(p1->c==p2->c){
	    a += SWF_MIN(p1->n,p2->n);
	    b += SWF_MAX(p1->n,p2->n);
	    p1 = p1->next;
	    p2 = p2->next;
	    continue;
	}
	if(p1->c < p2->c){
	    b += p1->n;
	    p1 = p1->next;
	    continue;
	}
	if(p1->c > p2->c){
	    b += p2->n;
	    p2 = p2->next;
	    continue;
	}
    }
    if(p1!=NULL){
	while(p1!=NULL){
	    b += p1->n;
	    p1 = p1->next;
	}
    }
    if(p2!=NULL){
	while(p2!=NULL){
	    b += p2->n;
	    p2 = p2->next;
	}
    }
    float ans;
    ans=3.0*a/b-1.0;
    return ans;
}

/*void *opt is a float pointer here, the value is the baseline*/
float MatchScore_Bow_Resemblance_ByNumber_Baseline(struct word_node *bow1, struct word_node *bow2, void *opt)
{
    struct word_node *p1;
    struct word_node *p2;
    p1 = bow1->next;
    p2 = bow2->next;
    /*a is number of intersection and b is union*/
    float bl;
    bl = *((float*)opt);
    float a,b;
    a=0.0;
    b=0.0;
    while(p1!=NULL&&p2!=NULL){
	if(p1->c==p2->c){
	    a += SWF_MIN(p1->n,p2->n) + bl;
	    b += SWF_MAX(p1->n,p2->n) + bl;
	    p1 = p1->next;
	    p2 = p2->next;
	    continue;
	}
	if(p1->c < p2->c){
	    b += p1->n + bl;
	    p1 = p1->next;
	    continue;
	}
	if(p1->c > p2->c){
	    b += p2->n + bl;
	    p2 = p2->next;
	    continue;
	}
    }
    if(p1!=NULL){
	while(p1!=NULL){
	    b += p1->n + bl;
	    p1 = p1->next;
	}
    }
    if(p2!=NULL){
	while(p2!=NULL){
	    b += p2->n + bl;
	    p2 = p2->next;
	}
    }
    float ans;
    ans=3.0*a/b-1.0;
    return ans;
}

void SWA_Bow_Even(struct word_node* bow_seq1, struct word_node* bow_seq2, int n1, int n2, int w, MatchingFunction_BOW MFB, float alpha, void* opt, float** map_score, unsigned char** map_trace)
{
    int i,j;
    int np1=n1+1;
    int np2=n2+1;
    int npw1=n1+w;
    int npw2=n2+w;
    int wminus1=w-1;
    /*Initialize*/
    for(i=0;i<npw1;i++){
	for(j=0;j<w;j++){
	    map_score[j][i]=0;
	    map_trace[j][i]=0;
	}
    }
    for(i=0;i<w;i++){
	for(j=0;j<npw2;j++){
	    map_score[j][i]=0;
	    map_trace[j][i]=0;
	}
    }
    /*update the matrix*/
    for(i=w;i<npw2;i++){
        /*printf("%d\n",i);*/
	int iminusw=i-w;
	int iminus1=i-1;
	for(j=w;j<npw1;j++){
	    int jminusw=j-w;
	    int jminus1=j-1;
	    float scoretemp;
	    float matchingscore;
	    matchingscore = MFB(bow_seq1+jminusw,bow_seq2+iminusw,opt);
            map_score[i][j] = 0;
	    map_trace[i][j] = 0;
	    scoretemp = matchingscore+map_score[iminusw][jminusw];
	    if(scoretemp>map_score[i][j]){
		map_score[i][j] = scoretemp;
		map_trace[i][j] = 1;
	    }
	    scoretemp = map_score[i][jminus1] - alpha;
	    if(scoretemp>map_score[i][j]){
		map_score[i][j] = scoretemp;
		map_trace[i][j] = 2;
	    }
	    scoretemp = map_score[iminus1][j] - alpha;
	    if(scoretemp>map_score[i][j]){
		map_score[i][j] = scoretemp;
		map_trace[i][j] = 3;
	    }
	}
    }
    return;
}


void Free_Bow(struct word_node *bow)
{
    struct word_node *node;
    struct word_node *last;
    last = bow->next;
    while(last!=NULL){
	node = last;
	last = node->next;
	free(node);
    }
    return;
}

void Bow_PlusNumber(struct word_node *bow_in, float addvalue, struct word_node *bow_out)
{
    if(bow_out==bow_in){
	struct word_node *nodethis;
	nodethis = bow_in->next;
	while(nodethis!=NULL){
	    nodethis->n += addvalue;
	}
    }
    else{
	struct word_node *nodethis;
	struct word_node *nodenew;
	struct word_node *nodelast;
	nodethis = bow_in->next;
	nodelast = bow_out;
	while(nodethis!=NULL){
	    nodenew = (struct word_node*)malloc(sizeof(struct word_node));
	    nodenew->n += nodethis->n+addvalue;
	    nodenew->c = nodethis->c;
	    nodenew->next = NULL;
	    nodelast->next = nodenew;
	    nodelast = nodenew;
	}
    }
    return;
}

/*Add 2 bow structure into one bow structure*/
void Bow_Plus(struct word_node *bow_in1, struct word_node *bow_in2, struct word_node *bow_out)
{
    struct word_node nodeout;
    struct word_node *bow_outtemp=NULL;
    if(bow_out==bow_in1||bow_out==bow_in2){
	bow_outtemp = bow_out;
	bow_out = &nodeout;
    }
    struct word_node *lasto;
    struct word_node *nodenew;
    struct word_node *p1;
    struct word_node *p2;
    lasto = bow_out;
    p1 = bow_in1->next;
    p2 = bow_in2->next;
    while(p1!=NULL&&p2!=NULL){
	if(p1->c==p2->c){
	    nodenew = (struct word_node*)malloc(sizeof(struct word_node));
	    nodenew->c = p1->c;
	    nodenew->n = p1->n+p2->n;
	    nodenew->next = NULL;
	    lasto->next = nodenew;
	    lasto = nodenew;
	    p1 = p1->next;
	    p2 = p2->next;
	    continue;
	}
	if(p1->c<p2->c){
	    nodenew = (struct word_node*)malloc(sizeof(struct word_node));
	    nodenew->c = p1->c;
	    nodenew->n = p1->n;
	    nodenew->next = NULL;
	    lasto->next = nodenew;
	    lasto = nodenew;
	    p1 = p1->next;
	    continue;
	}
	if(p1->c>p2->c){
	    nodenew = (struct word_node*)malloc(sizeof(struct word_node));
	    nodenew->c = p2->c;
	    nodenew->n = p2->n;
	    nodenew->next = NULL;
	    lasto->next = nodenew;
	    lasto = nodenew;
	    p2 = p2->next;
	    continue;
	}
    }
    if(p1!=NULL){
	while(p1!=NULL){
	    nodenew = (struct word_node*)malloc(sizeof(struct word_node));
	    nodenew->c = p1->c;
	    nodenew->n = p1->n;
	    nodenew->next = NULL;
	    lasto->next = nodenew;
	    lasto = nodenew;
	    p1 = p1->next;
	}
    }
    if(p2!=NULL){
	while(p2!=NULL){
	    nodenew = (struct word_node*)malloc(sizeof(struct word_node));
	    nodenew->c = p2->c;
	    nodenew->n = p2->n;
	    nodenew->next = NULL;
	    lasto->next = nodenew;
	    lasto = nodenew;
	    p2 = p2->next;
	}
    }
    if(bow_outtemp==NULL){
	return;
    }
    else{
        if(bow_outtemp==bow_in1){
	    Free_Bow(bow_in1);
	    bow_in1->next = nodeout.next;
	    return;
	}
	else{
	    Free_Bow(bow_in2);
	    bow_in2->next = nodeout.next;
	    return;
	}
    }
}

/*bow3=bow1-bow2*/
void Bow_Minus(struct word_node* bow_in1, struct word_node* bow_in2, struct word_node* bow_out)
{
    struct word_node nodeout;
    struct word_node *bow_outtemp=NULL;
    if(bow_out==bow_in1||bow_out==bow_in2){
	bow_outtemp = bow_out;
	bow_out = &nodeout;
    }
    struct word_node *lasto;
    struct word_node *nodenew;
    struct word_node *p1;
    struct word_node *p2;
    lasto = bow_out;
    p1 = bow_in1->next;
    p2 = bow_in2->next;
    while(p1!=NULL&&p2!=NULL){
	if(p1->c==p2->c){
	    if( p1->n - p2->n > Fre_Zero ){
		nodenew = (struct word_node*)malloc(sizeof(struct word_node));
		nodenew->c = p1->c;
		nodenew->n = p1->n-p2->n;
		nodenew->next = NULL;
		lasto->next = nodenew;
		lasto = nodenew;
	    }
	    p1 = p1->next;
	    p2 = p2->next;
	    continue;
	}
	if(p1->c<p2->c){
	    nodenew = (struct word_node*)malloc(sizeof(struct word_node));
	    nodenew->c = p1->c;
	    nodenew->n = p1->n;
	    nodenew->next = NULL;
	    lasto->next = nodenew;
	    lasto = nodenew;
	    p1 = p1->next;
	    continue;
	}
	if(p1->c>p2->c){
	    p2 = p2->next;
	    continue;
	}
    }
    if(p1!=NULL){
	while(p1!=NULL){
	    nodenew = (struct word_node*)malloc(sizeof(struct word_node));
	    nodenew->c = p1->c;
	    nodenew->n = p1->n;
	    nodenew->next = NULL;
	    lasto->next = nodenew;
	    lasto = nodenew;
	    p1 = p1->next;
	}
    }
    if(bow_outtemp==NULL){
	return;
    }
    else{
        if(bow_outtemp==bow_in1){
	    Free_Bow(bow_in1);
	    bow_in1->next = nodeout.next;
	    return;
	}
	else{
	    Free_Bow(bow_in2);
	    bow_in2->next = nodeout.next;
	    return;
	}
    }
}

/*Alignment pairs are malloced inside this function*/
/* *align is father node, malloced outside.*/
void Trace_Bow_Even(unsigned char **map_trace, int w, int p1, int p2, struct pair_node *align)
{
    align->p1 = p1;
    align->p2 = p2;
    align->next = NULL;
    struct pair_node *pn;
    while(1){
	if(map_trace[p2][p1]==0)
	    break;
	if(map_trace[p2][p1]==1){
	    pn = (struct pair_node*)malloc(sizeof(struct pair_node));
	    p2 -= w;
	    p1 -= w;
	    pn->p1 = p1;
	    pn->p2 = p2;
	    pn->next = align->next;
	    align->next = pn;
	    continue;
	}
	if(map_trace[p2][p1]==2){
	    p1 -= 1;
	    continue;
	}
	if(map_trace[p2][p1]==3){
	    p2 -= 1;
	    continue;
	}
    }
    return;
}


void Print_Alignment_Bow(struct pair_node *align, unsigned char *seq1, unsigned char *seq2, int m, int k, FILE *outfile)
{
    if(align->next==NULL){
	return;
    }
    int w=m/k;
    int p1,p2;
    p1 = align->next->p1*k;
    p2 = align->next->p2*k;
    int s1,t1;
    int s2,t2;
    struct pair_node *pn;
    pn = align->next;
    int i;
    fprintf(outfile,"\n%d",p2);
    while(pn!=NULL){
	s1 = pn->p1*k;
	s2 = pn->p2*k;
	for(i=p1;i<s1;i++){
	    fprintf(outfile,"\n%c,%c",State2Char(seq1[i]),'-');
	}
	for(i=p2;i<s2;i++){
	    fprintf(outfile,"\n%c,%c",'-',State2Char(seq2[i]));
	}
        fprintf(outfile,"\n(");
	for(i=0;i<m;i++){
	    fprintf(outfile,"%c ",State2Char(seq1[s1+i]));
	}
	fprintf(outfile,")  (");
	for(i=0;i<m;i++){
	    fprintf(outfile,"%c ",State2Char(seq2[s2+i]));
	}
	fprintf(outfile,")");
	p1 = s1+m;
	p2 = s2+m;
        pn = pn->next;
    }
    fprintf(outfile,"\n");
    return;
}


void Trace_Even(unsigned char **map_trace, int p1, int p2, struct pair_node *align)
{
    align->p1 = p1;
    align->p2 = p2;
    align->next = NULL;
    struct pair_node *pn;
    while(1){
	if(map_trace[p2][p1]==0)
	    break;
	if(map_trace[p2][p1]==1){
	    pn = (struct pair_node*)malloc(sizeof(struct pair_node));
	    p1 -= 1;
	    p2 -= 1;
	    pn->p1 = p1;
	    pn->p2 = p2;
	    pn->next = align->next;
	    align->next = pn;
	    continue;
	}
	if(map_trace[p2][p1]==2){
	    p1 -= 1;
	    continue;
	}
	if(map_trace[p2][p1]==3){
	    p2 -= 1;
	    continue;
	}
    }
    return;
}

void Print_Alignment_Even(struct pair_node *align, unsigned char *seq1, unsigned char *seq2, FILE *outfile)
{
    int i;
    int p1,p2;
    int s1,s2;
    align = align->next;
    p1 = align->p1;
    p2 = align->p2;
    while(align!=NULL){
	s1 = align->p1;
	s2 = align->p2;
	for(i=p1;i<s1;i++){
	    fprintf(outfile,"%c,%c\n",State2Char(seq1[i]),'-');
	}
	for(i=p2;i<s2;i++){
	    fprintf(outfile,"%c,%c\n",'-',State2Char(seq2[i]));
	}
	fprintf(outfile,"%c,%c\n",State2Char(seq1[s1]),State2Char(seq2[s2]));
        p1 = s1+1;
        p2 = s2+1;
	align = align->next;
    }
    return;
}

void Print_Alignment_Sseq_Even(struct pair_node *align, unsigned char *sseq1, unsigned char *sseq2, unsigned short *sseq1_num, unsigned short *sseq2_num, FILE *outfile)
{
    int i;
    int p1,p2;
    int s1,s2;
    align = align->next;
    p1 = align->p1;
    p2 = align->p2;
    while(align!=NULL){
	s1 = align->p1;
	s2 = align->p2;
	for(i=p1;i<s1;i++){
	    fprintf(outfile,"%c,%c\t%d,%d\n",State2Char(sseq1[i]),'-',sseq1_num[i],0);
	}
	for(i=p2;i<s2;i++){
	    fprintf(outfile,"%c,%c\t%d,%d\n",'-',State2Char(sseq2[i]),0,sseq2_num[i]);
	}
	fprintf(outfile,"%c,%c\t%d,%d\n",State2Char(sseq1[s1]),State2Char(sseq2[s2]),sseq1_num[s1],sseq2_num[s2]);
        p1 = s1+1;
        p2 = s2+1;
	align = align->next;
    }
    return;
}

void Free_Alignment(struct pair_node *align)
{
    struct pair_node *pn;
    align = align->next;
    while(align!=NULL){
        pn = align;
	align = align->next;
	free(pn);
    }
    return;
}

/*layer 4(value of map_trace_z) is the empty layer*/
void SWA_Linear(unsigned char *seq1, unsigned char *seq2, int l1, int l2, MatchingFunction MF, GapFunction GF, float alpha, float beta, void *opt, float ***map_score, unsigned char ***map_trace)
{
    int lp1;
    int lp2;
    lp1 = l1+1;
    lp2 = l2+1;
    int i,j,k;
    /*Initialize*/
    for(k=0;k<4;k++){
	for(i=0;i<lp1;i++){
	    map_score[k][0][i] = 0.0;
	    map_trace[k][0][i] = 4;
	}
    }
    for(k=0;k<4;k++){
	for(j=0;j<lp2;j++){
	    map_score[k][j][0] = 0.0;
	    map_trace[k][j][0] = 4;
	}
    }
    /*Do dynamic programming*/
    int zt;
    float score;
    int ip,jp;
    for(j=0;j<lp2;j++){
	jp = j+1;
	for(i=0;i<lp1;i++){
	    ip = i+1;
	    /*k = 0*/
	    score = MF(seq1[i],seq2[j],i,j,opt);
	    /*k = 1*/
	    if(map_score[0][jp][i]>map_score[1][jp][i]){
	    }
	    else{
	    }
	    /*k = 2*/
            if(map_score[0][j][ip]>map_score[2][j][ip]){
	    }
	    else{
	    }
	    /*k = 3*/
	}
    }
    return;
}

void Trace_Linear(unsigned char ***map_trace, int p1, int p2, int z, struct pair_node *align)
{
}



