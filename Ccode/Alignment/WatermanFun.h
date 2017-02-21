#ifndef WatermanFun_H
#define WatermanFun_H

#include<stdlib.h>
#include<stdio.h>
#include<limits.h>
#include"EpiBLASTConstant.h"

#define SWF_MAX(a, b) (( a > b) ? a : b)
#define SWF_MIN(a, b) (( a > b) ? b : a)

struct Map_State_Even
{
    int s;
    int l1;
    float *score;
};

/*To define bag of words structure*/
struct word_node
{
    unsigned char c;
    float n;
    struct word_node *next;
};

/*Facilitate bag of words computation from raw sequences*/
struct word_node_extend
{
    unsigned char c;
    float n;
    int p;
    struct word_node_extend *next;
};

struct pair_node
{
    int p1;
    int p2;
    struct pair_node *next;
};

typedef float (*MatchingFunction)(unsigned char, unsigned char, int, int, void*);
typedef float (*GapFunction)(float, int, char, void*);
typedef float (*MatchingFunction_BOW)(struct word_node*, struct word_node*, void*);

/*Convert chromatin state to human comprehensible characters*/
char State2Char(unsigned char);
/*Reverse function of State2Char*/
unsigned char Char2State(char);
/*Convert state sequence to compact sequence, output arrays are malloced inside this function*/
int Seq2Sseq(unsigned char*, int, unsigned char**, unsigned short**);
/*Convert compact sequence to state sequence, output array is malloced inside this function*/
int Sseq2Seq(unsigned char*, unsigned short*, int, unsigned char **);

float MatchScore_Naive(unsigned char, unsigned char, int, int, void*);
float MatchScore_Sqrt(unsigned char, unsigned char, int, int, void*);
float GapScore_Naive(float, int, char, void*);
/*In this function, gap penalty is alpha*l*/
void SWA_Even(unsigned char*, unsigned char*, int, int, MatchingFunction, GapFunction, float, struct Map_State_Even*, void*, float**, unsigned char**);
void Trace_Even(unsigned char**, int, int, struct pair_node*);
void Print_Alignment_Even(struct pair_node*, unsigned char*, unsigned char*, FILE *);
void Print_Alignment_Sseq_Even(struct pair_node*, unsigned char*, unsigned char*, unsigned short*, unsigned short*, FILE *);

/*In this function, gap penalty is alpha*l+beta*/
void SWA_Linear(unsigned char*, unsigned char*, int, int, MatchingFunction, GapFunction, float, float, void*, float***, unsigned char***);
void Trace_Linear(unsigned char***, int, int, int, struct pair_node*);

void Malloc_Map_Compact(unsigned char*, unsigned char*, unsigned short*, unsigned short*, int, int, float***, unsigned char ***);
void Free_Map(int, float**, unsigned char **);
void Free_Alignment(struct pair_node*);

void SWA_Compact_Even(unsigned char*, unsigned char*, unsigned short*, unsigned short*, int, int, MatchingFunction, float, void*, float**, unsigned char**);


/*bow_seq is malloced inside this function*/
/*state 0 is reserved state here*/
/*Seq2Bow_Slow does the same thing with function Seq2Bow*/
/*Seq2Bow_Slow is slow when m/k was small, but Seq2Bow_Slow is insensitive to m/k hence have better speed when m/k was large*/
int Seq2Bow(unsigned char *, int, int, int, struct word_node**);
int Seq2Bow_Slow(unsigned char *, int, int, int, struct word_node**);

/*These functions define basic operations of bag of words structure*/
/*If the output pointer doesn't equal to the input pointer, new space is malloced, otherwise the corresponding input structure is destroyed and re-assigned*/
/*Add 2 bow structure into one bow structure*/
void Bow_Plus(struct word_node*, struct word_node*, struct word_node*);
/*bow3=bow1-bow2*/
void Bow_Minus(struct word_node*, struct word_node*, struct word_node*);
/*A number added to all the words exist in a bow structure*/
void Bow_PlusNumber(struct word_node*, float, struct word_node*);

/*This function releases the space of a bag of words structure*/
void Free_Bow(struct word_node*);

float MatchScore_Bow_Resemblance(struct word_node*, struct word_node*, void*);
float MatchScore_Bow_Resemblance_ByNumber(struct word_node*, struct word_node*, void*);
/*void *opt is a float pointer here, the value is the baseline*/
float MatchScore_Bow_Resemblance_ByNumber_Baseline(struct word_node*, struct word_node*, void*);

void SWA_Bow_Even(struct word_node*, struct word_node*, int, int, int, MatchingFunction_BOW, float, void*, float**, unsigned char**);

/*Alignment pairs are malloced inside this function*/
/* *align is father node, malloced outside.*/
void Trace_Bow_Even(unsigned char **, int, int, int, struct pair_node*);

void Print_Alignment_Bow(struct pair_node*, unsigned char*, unsigned char*, int, int, FILE *);

#endif




