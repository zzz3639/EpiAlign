#include<stdio.h>
#include<time.h>
#include"WatermanFun.h"

void readfile(int n, const char *filename, unsigned char *sseq, int *sseq_num)
{
    int i;
    FILE *infile;
    int a,b;
    infile=fopen(filename,"r");
    for(i=0;i<n;i++){
	fscanf(infile,"%d",&a);
	fscanf(infile,"%d",&b);
	sseq[i]=a;
	sseq_num[i]=b;
    }
    fclose(infile);
    return;
}

int main(int argc, char **argv)
{
    int n;
    n = atoi(argv[2]);
    unsigned char *sseq;
    int *sseq_num;
    sseq = (unsigned char *)malloc(sizeof(unsigned char)*n);
    sseq_num = (int *)malloc(sizeof(int)*n);
    readfile(n,argv[1],sseq,sseq_num);

    int m=0;
    int i,j;
    for(i=0;i<n;i++){
	m+=sseq_num[i];
    }

    unsigned char *seq;
    seq = (unsigned char *)malloc(sizeof(unsigned char)*m);
    m=0;
    for(i=0;i<n;i++){
	for(j=0;j<sseq_num[i];j++){
	    seq[m]=sseq[i];
	    m+=1;
	}
    }
    n=m;

    struct word_node *bow_seq1;
    struct word_node *bow_seq2;
    m = atoi(argv[3]);
    int k=1;
    int l1;
    int l2;
    clock_t t1,t2,t3;
    t1 = clock();
    l1 = Seq2Bow_Slow(seq, n, m, k, &bow_seq1);
    t2 = clock();
    l2 = Seq2Bow(seq, n, m, k, &bow_seq2);
    t3 = clock();
    printf("\n%d, %d\n",t2-t1,t3-t2);

    int eq=1;
    if(l1!=l2){
	printf("\nNot equal!\n");
	return 1;
    }
    struct word_node *wordnode1;
    struct word_node *wordnode2;
    for(i=0;i<l1;i++){
	wordnode1 = bow_seq1[i].next;
	wordnode2 = bow_seq2[i].next;
	while(wordnode1!=NULL&&wordnode2!=NULL){
	    if(wordnode1->c!=wordnode2->c){
		printf("\nNot equal!\n");
		return 1;
	    }
	    if(wordnode1->n!=wordnode2->n){
		printf("\nNot equal!\n");
		return 1;
	    }
	    wordnode1 = wordnode1->next;
	    wordnode2 = wordnode2->next;
	}
	if(wordnode1!=NULL||wordnode2!=NULL){
	    printf("\nNot equal!\n");
	    return 1;
	}
    }

    printf("\nequal\n");
    return 1;
}




