#ifndef StateIO_H
#define StateIO_H

#include<stdio.h>
#include<stdlib.h>
#include<ctype.h>
#include<limits.h>
#include"EpiBLASTConstant.h"

#define ASCII_MAX 128
#define FileNameLength_MAX 200
#define LineStore_MAX 200

struct StateIO_Opt
{
    char s;
    int n;
    char key[ASCII_MAX];
};

struct Matrix_Size
{
    /*line number*/
    int m;
    /*column number*/
    int n;
};

/*please be noticed, struct StateIO_Opt *opt should be malloced outside of this function*/
void MakeKey_Number(struct StateIO_Opt *);
void MakeKey_Separator(struct StateIO_Opt *);

/*For struct StateIO_Opt *opt*/
/*If opt->s==0 or opt==NULL, then lines are counted*/
/*If opt->s==1, then lines started as defined by opt->key are counted*/
int FILE_CountLine(const char*, struct StateIO_Opt *);

/*For struct StateIO_Opt *opt*/
/*If opt==NULL, then opt is assigned inside this function with opt->s=0 and opt->n<=0*/
/*opt->key should have been assigned for number charactors*/
/*If opt->s==0, then Sseq and Sseq_num are malloced inside*/
/*If opt->s==1, then Sseq and Sseq_num are malloced outside*/
/*If opt->n<=0, then line number is re-counted*/
/*If opt->n>0, then line number is set to opt->n*/
int Sseq_ReadFile(const char*, unsigned char**, unsigned short**, struct StateIO_Opt *);

/*For struct StateIO_Opt *opt*/
/*If opt==NULL, then opt is assigned inside this function with opt->s=0 and opt->n<=0*/
/*opt->key should have been assigned for number charactors*/
/*If opt->s==0, then Sseq and Sseq_num are malloced inside*/
/*If opt->s==1, then Sseq and Sseq_num are malloced outside*/
/*If opt->n<=0, then line number is re-counted*/
/*If opt->n>0, then line number is set to opt->n*/
int Seq_ReadFile(const char*, unsigned char**, struct StateIO_Opt *);

/*ToBeFinished*/
void Matrix_Parsing(const char*, struct StateIO_Opt *, struct Matrix_Size *);

/*ToBeFinished*/
void Posterior_ReadFile(const char*, float***, struct StateIO_Opt *);

/*ToBeFinished*/
/*If opt==NULL, then lines are counted inside, lines container is malloced inside*/
int Lines_ReadFile(const char*, char ***, struct StateIO_Opt *);

#endif

