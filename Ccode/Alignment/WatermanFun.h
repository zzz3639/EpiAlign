#include<stdlib.h>
#include<limits.h>

#define STATE_MAX 256
#define SWF_MAX(a, b) (( a > b) ? a : b)
#define SWF_MIN(a, b) (( a > b) ? a : b)

struct Map_State_Even
{
    int s;
    int l1;
    float *score;
};

struct word_node
{
    unsigned char c;
    unsigned short n;
    struct word_node *next;
};

struct word_node_extend
{
    unsigned char c;
    unsigned short n;
    int p;
    struct word_node_extend *next;
};

typedef float (*MatchingFunction)(unsigned char, unsigned char, void*);
typedef float (*MatchingFunction_BOW)(struct word_node*, struct word_node*, void*);

char State2Char(unsigned char);
/*Convert chromatin state to human comprehensible characters*/
unsigned char Char2State(char);
/*Reverse function of State2Char*/
int Seq2Sseq(unsigned char*, int, unsigned char**, unsigned short**);
/*Convert state sequence to compact sequence, output arrays are malloced inside this function*/
float MatchScore_Naive(unsigned char, unsigned char, void*);
void SWA_Even(unsigned char*, unsigned char*, int, int, MatchingFunction, float, struct Map_State_Even*, void*, float**, unsigned char**);
void Trace_Even(void);
void Malloc_Map_Compact(unsigned char*, unsigned char*, unsigned short*, unsigned short*, int, int, float***, unsigned char ***);
void Free_Map(int, float**, unsigned char **);
void SWA_Compact_Even(unsigned char*, unsigned char*, unsigned short*, unsigned short*, int, int, MatchingFunction, float, void*, float**, unsigned char**);


/*bow_seq is malloced inside this function*/
/*state 0 is reserved state here*/
/*Seq2Bow_Slow does the same thing with function Seq2Bow, but time cost of this function is higher*/
int Seq2Bow_Slow(unsigned char *, int, int, int, struct word_node**);
int Seq2Bow(unsigned char *, int, int, int, struct word_node**);

void Bow_Plus(struct word_node*, struct word_node*);
void Bow_Minus(struct word_node*, struct word_node*);
void Free_Bow(struct word_node*);

float MatchScore_Bow_Resemblance(struct word_node*, struct word_node*, void*);
float MatchScore_Bow_Resemblance_ByNumber(struct word_node*, struct word_node*, void*);
void SWA_Bow_Even(struct word_node*, struct word_node*, int, int, int, MatchingFunction_BOW, float, void*, float**, unsigned char**);




