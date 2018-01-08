#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"StateIO.h"

#define chrnamekeylen 3
#define chrlettermax 3

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

void QuickSort(int *A, int *idx, int m, int n)
{
    int key;
    int key_idx;
    int i,j,k;
    if(m<n){
	k = ChoosePivot(m,n);
	SwapInt(A+m,A+k);
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
		SwapInt(A+i,A+j);
		SwapInt(idx+i,idx+j);
	    }
	}
	SwapInt(A+m,A+j);
	SwapInt(idx+m,idx+j);
	QuickSort(A,idx,m,j-1);
	QuickSort(A,idx,j+1,n);
    }
    return;
}

int ChrStr2num(char *str, int N)
{
    int ans;
    if(str[0]>='0'&&str[0]<='9'){
	ans = atoi(str);
    }
    else{
	ans = str[0]+N;
    }
    return ans;
}

int main(int argc, char **argv)
{
    if(argc==1){
	printf("\nUsage: run.out indexin indexout\n");
	return 1;
    }
/*read input*/
    int l;
    char **chrstrings;
    l = Lines_ReadFile(argv[1], &chrstrings, NULL);
/*find the identity string for each chromosome name*/
    char chrnamekey[chrnamekeylen+1];
    chrnamekey[0]='c'; chrnamekey[1]='h'; chrnamekey[2]='r'; chrnamekey[3]='\0';
    char **chrletters;
    chrletters = (char**)malloc(sizeof(char*)*l);
    int i,j;
    for(i=0;i<l;i++){
	chrletters[i] = (char*)malloc(sizeof(char)*(chrlettermax+1));
    }
    char *strsearchtemp;
    for(i=0;i<l;i++){
	strsearchtemp = strstr(chrstrings[i],chrnamekey);
	j=chrnamekeylen;
	while(strsearchtemp[j]!='_'){
	    chrletters[i][j-chrnamekeylen] = strsearchtemp[j];
	    j++;
	}
	chrletters[i][j-chrnamekeylen] = '\0';
    }
/*sort the identity strings*/
    int *chrvalues;
    int *chrindexes;
    chrvalues = (int*)malloc(sizeof(int)*l);
    chrindexes = (int*)malloc(sizeof(int)*l);
    for(i=0;i<l;i++){
	chrindexes[i] = i;
	chrvalues[i] = ChrStr2num(chrletters[i], l);
    }
    QuickSort(chrvalues,chrindexes,0,l-1);
/*write to the output*/
    FILE *fout;
    fout = fopen(argv[2],"w");
    for(i=0;i<l;i++){
	fprintf(fout,"%s\n",chrstrings[chrindexes[i]]);
    }
    fclose(fout);
/*free malloced spaces*/
    for(i=0;i<l;i++){
	free(chrstrings[i]);
	free(chrletters[i]);
    }
    free(chrstrings);
    free(chrletters);
    free(chrvalues);
    free(chrindexes);
    return 1;
}



