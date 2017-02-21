/*This code generates fake chromosome based on real chromosome by 1-gram*/
/*Compile this code by g++4.8 or later version, as this code uses c++11 grammer*/
#include<stdio.h>
#include<stdlib.h>
#include"StateIO.h"
#include"WatermanFun.h"
#include<random>
#include<vector>
#include<chrono>

int main(int argc, char **argv)
{
    if(argc!=3){
	printf("\nUsage: ./run chromosome fake_chromosome\n");
	return 1;
    }
/*read true chromosome*/
    unsigned char *sseq;
    unsigned short *sseq_num;
    int nt;
    nt = Sseq_ReadFile(argv[1],&sseq,&sseq_num,NULL);
    int lt = 0;
    int i,j,k;
    for(i=0;i<nt;i++){
	lt += sseq_num[i];
    }
    if(lt==0){
	FILE *outf = fopen(argv[2],"w");
	fclose(outf);
	return 1;
    }
    if(lt==1){
	FILE *outf = fopen(argv[2],"w");
	fprintf(outf,"%d %d\n",sseq[0],sseq_num[0]);
	fclose(outf);
	return 1;
    }
    unsigned char *seq;
    seq = (unsigned char*)malloc(sizeof(unsigned char)*lt);
    k=0;
    for(i=0;i<nt;i++){
	for(j=0;j<sseq_num[i];j++){
	    seq[k] = sseq[i];
	    k++;
	}
    }
    char *StateCount;
    StateCount = (char*)malloc(sizeof(char)*STATE_MAX);
    for(i=0;i<STATE_MAX;i++){
	StateCount[i] = 0;
    }
    for(i=0;i<nt;i++){
        StateCount[sseq[i]] = 1;
    }
    int ns=0;
    for(i=0;i<STATE_MAX;i++){
	if(StateCount[i]){
	    ns = i+1;
	}
    }
/*compute frequency and 1-gram frequency*/
    std::vector<double> fre(ns,0.0);
    std::vector<std::vector<double> > g1_fre(ns,std::vector<double>(ns,0.0));
    for(i=0;i<lt;i++){
	fre[seq[i]] += 1.0;
    }
    for(i=0;i<lt-1;i++){
	g1_fre[seq[i]][seq[i+1]] += 1.0;
    }
    /*printf("\n");
    for(i=0;i<ns;i++){
        printf("%f ",fre[i]);
    }
    printf("\n\n");*/
    /*for(i=0;i<ns;i++){
	for(j=0;j<ns;j++){
	    printf("%f ",g1_fre[i][j]);
	}
	printf("\n");
    }*/
/*Initialize random number generator*/
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937_64 gen(seed);
/*generate fake chromosome*/
    unsigned char *seq_fake;
    seq_fake = (unsigned char *)malloc(sizeof(unsigned char)*lt);
    std::discrete_distribution<int> MN_Gen_first(fre.begin(),fre.end());
    seq_fake[0] = MN_Gen_first(gen);
    std::vector<std::discrete_distribution<int> > MN_Gen_g1(ns,std::discrete_distribution<int>());
    for(i=0;i<ns;i++){
	MN_Gen_g1[i] = std::discrete_distribution<int>(g1_fre[i].begin(),g1_fre[i].end());
    }
    for(i=1;i<lt;i++){
	seq_fake[i] = MN_Gen_g1[seq_fake[i-1]](gen);
    }
    unsigned char *sseq_fake;
    unsigned short *sseq_fake_num;
    int no;
    no = Seq2Sseq(seq_fake,lt,&sseq_fake,&sseq_fake_num);

/*output result*/
    FILE *outf = fopen(argv[2],"w");
    for(i=0;i<no;i++){
	fprintf(outf,"%d %d\n",sseq_fake[i],sseq_fake_num[i]);
    }
    fclose(outf);
/*free memory*/
    free(seq);
    free(sseq);
    free(sseq_num);
    free(StateCount);
    free(seq_fake);
    free(sseq_fake);
    free(sseq_fake_num);
    return 1;
}


