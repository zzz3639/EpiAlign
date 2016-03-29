#include"WatermanFun.h"

char State2Char(unsigned char c)
{
    return (c<27)?(char)(c)+'`':(char)(c)+'&';
}

unsigned char Char2State(char c)
{
    return (c>'`')?(unsigned char)c-(unsigned char)96:(unsigned char)c-(unsigned char)38;
}

float MatchScore_Naive(unsigned char a, unsigned char b, void* opt)
{
    return (a==b)?2.0:-1.0;
}

/*Smith-Waterman function for even gap penalty*/
/*By default, seq1 is query and seq2 is database*/
void SWA_Even(unsigned char* seq1, unsigned char* seq2, int l1, int l2, MatchingFunction MF, float alpha, struct Map_State_Even* map_init, void* opt, float** map_score, unsigned char** map_trace)
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
	    float matchthis=MF(seq1[jloop],seq2[iloop],opt);
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
	    scoretemp = map_score[i][jloop] - alpha;
	    if(scoretemp>map_score[i][j]){
		map_score[i][j] = scoretemp;
		map_trace[i][j] = 2;
	    }
	    /*deletion in sequence 2*/
	    scoretemp = map_score[iloop][j] - alpha;
	    if(scoretemp>map_score[i][j]){
		map_score[i][j] = scoretemp;
		map_trace[i][j] = 3;
	    }
	}
    }
    return;
}

void Seq2Sseq(unsigned char *seq, int n, unsigned char *sseq, unsigned short *sseq_len)
{
    if(n<1){
	sseq=NULL;
	sseq_len=NULL;
	return;
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
    sseq=(unsigned char *)malloc(sizeof(unsigned char)*k);
    sseq_len=(unsigned short *)malloc(sizeof(unsigned short)*k);
    s=seq[0];
    t=1;
    j=0;
    for(i=1;i<n;i++){
	if(seq[i]==s){
	    t+=1;
	    continue;
	}
	if(t<=USHRT_MAX){
	    sseq[j]=s;
	    sseq_len[j]=t;
	    j+=1;
	}
	else{
	    r=(t+USHRT_m1)/USHRT_MAX;
            for(m=0;m<r-1;m++){
		sseq[j]=s;
		sseq_len[j]=USHRT_MAX;
		j+=1;
	    }
	    sseq[j]=s;
	    sseq_len[j]=(t-1)%USHRT_MAX+1;
	    j+=1;
	}
    }
    return;
}

void Malloc_Map_Compact(unsigned char* seq1, unsigned char* seq2, unsigned short* len1, unsigned short* len2, int l1, int l2, float** map_score, unsigned char** map_trace)
{
    int i,j,k;
    int lp1=l1+1;
    int lp2=l2+1;
    int n1=0,n2=0;
    for(i=0;i<l1;i++)
	n1+=len1[i];
    for(i=0;i<l2;i++)
	n2+=len2[i];
    map_score=(float**)malloc(sizeof(float*)*(n2+1));
    map_trace=(unsigned char**)malloc(sizeof(unsigned char*)*(n2+1));

    /*malloc the columns*/
    int np1=n1+1;
    char *M=(char*)malloc(sizeof(char)*(n2+1));
    for(i=0;i<=n2;i++)
	M[i]=0;
    map_score[0]=(float*)malloc(sizeof(float)*np1);
    map_trace[0]=(unsigned char*)malloc(sizeof(unsigned char)*np1);
    M[0]=1;
    k=0;
    for(i=0;i<l2;i++){
	k+=len2[i];
	map_score[k]=(float*)malloc(sizeof(float)*np1);
	map_trace[k]=(unsigned char *)malloc(sizeof(unsigned char)*np1);
        M[k]=1;
    }
    for(i=0;i<=n2;i++){
	if(M[i]==1)
	    continue;
	map_score[i]=(float*)malloc(sizeof(float)*lp1);
	map_trace[i]=(unsigned char*)malloc(sizeof(unsigned char)*lp1);
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
    r2=len2[0];
    t2=0;
    for(i=0;i<n2;i++){
	if(M[i]==0){
	    /*processing colume i item 0--(l1-1)*/
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
	    /*Processing colume i item l1*/
	    score_d2=map_score[i][l1] - r2*alpha;
	    if(score_d2>map_score[iplusr2][n1]){
		map_score[iplusr2][n1]=score_d2;
		map_trace[iplusr2][n1]=3;
	    }
	}
	else{
	    t1=0;
	    r1=len1[t1];
	    /*processing colume i item 0--(n1-1)*/
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
	    /*processing colume i item n1*/
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
    return;
}



