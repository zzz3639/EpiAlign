#include<stdio.h>
#include<stdlib.h>

#define NOTK_MAX(a, b) (( a > b) ? a : b)
#define NOTK_MIN(a, b) (( a > b) ? b : a)
#define lineofscores 2
#define MinusInf -1000000.0
#define MinusInfInt -1000000

void Remove_Peak(float **maxline, int *L, int h, int p, int dist_peak_l, int dist_peak_u)
{
    int l;
    int u;
    l = NOTK_MAX(p-dist_peak_l,0);
    u = NOTK_MIN(p+dist_peak_u,L[h]-1);
    int i;
    for(i=l;i<=u;i++){
	maxline[h][i] = MinusInf;
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

int IntMax_Array(int *data, int n)
{
    int ans=MinusInfInt;
    int i;
    for(i=0;i<n;i++){
	if(data[i]>ans)
	    ans = data[i];
    }
    return ans;
}

int main(int argc, char **argv)
{
    if(argc==1){
	printf("\nUsage: run.out #lines window K <scores >outputfile\n");
	return 1;
    }
    int M;
    int W;
    int K;
    M = atoi(argv[1]);
    W = atoi(argv[2]);
    K = atoi(argv[3]);
    float **scores;
    int *chr_full;
    int *pos_full;
    scores = (float**)malloc(sizeof(float*)*lineofscores);
    chr_full = (int*)malloc(sizeof(int)*M);
    pos_full = (int*)malloc(sizeof(int)*M);
    int i,j;
    for(i=0;i<lineofscores;i++){
	scores[i] = (float*)malloc(sizeof(float)*M);
    }
    for(i=0;i<M;i++){
	fscanf(stdin,"%f",scores[0]+i);
	fscanf(stdin,"%f",scores[1]+i);
	fscanf(stdin,"%d",chr_full+i);
	fscanf(stdin,"%d",pos_full+i);
    }
    int hmax;
    hmax = IntMax_Array(chr_full,M);
    hmax++;
    int *Lmatrix;
    Lmatrix = (int*)malloc(sizeof(int)*hmax);
    for(i=0;i<hmax;i++){
	Lmatrix[i] = 0;
    }
    for(i=0;i<M;i++){
	if(Lmatrix[chr_full[i]]<pos_full[i]){
	    Lmatrix[chr_full[i]] = pos_full[i];
	}
    }
    float **score_matrix;
    float **score_matrix_ori;
    score_matrix = (float**)malloc(sizeof(float*)*hmax);
    score_matrix_ori = (float**)malloc(sizeof(float*)*hmax);
    for(i=0;i<hmax;i++){
	score_matrix[i] = (float*)malloc(sizeof(float)*Lmatrix[i]);
	score_matrix_ori[i] = (float*)malloc(sizeof(float)*Lmatrix[i]);
    }
    for(i=0;i<M;i++){
	score_matrix[chr_full[i]][pos_full[i]-1] = (scores[1][i]-scores[0][i])/scores[0][i];
	score_matrix_ori[chr_full[i]][pos_full[i]-1] = score_matrix[chr_full[i]][pos_full[i]-1];
    }
    int chr_max;
    int pos_max;
    for(i=0;i<K;i++){
	Find_Peak(score_matrix, hmax, Lmatrix, &chr_max, &pos_max);
	Remove_Peak(score_matrix, Lmatrix, chr_max, pos_max, W-1, W-1);
	printf("%d\t%d\t%f\n",chr_max,pos_max+1,score_matrix_ori[chr_max][pos_max]);
    }
    return 1;
}






