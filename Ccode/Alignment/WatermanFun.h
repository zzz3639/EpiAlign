#include<stdlib.h>
#include<limits.h>

typedef float (*MatchingFunction)(unsigned char, unsigned char, void*);
struct Map_State_Even
{
    int s;
    int l1;
    float *score;
};

char State2Char(unsigned char);
/*Convert chromatin state to human comprehensible characters*/
unsigned char Char2State(char);
/*Reverse function of State2Char*/
void Seq2Sseq(unsigned char*, int, unsigned char*, unsigned short*);
/*Convert state sequence to compact sequence, output arrays are malloced inside this function*/
float MatchScore_Naive(unsigned char, unsigned char, void*);
void SWA_Even(unsigned char*, unsigned char*, int, int, MatchingFunction, float, struct Map_State_Even*, void*, float**, unsigned char**);
void Trace_Even(void);
void Malloc_Map_Compact(unsigned char*, unsigned char*, unsigned short*, unsigned short*, int, int, float**, unsigned char **);
void Free_Map(int, float**, unsigned char **);
void SWA_Compact_Even(unsigned char*, unsigned char*, unsigned short*, unsigned short*, int, int, MatchingFunction, float, void*, float**, unsigned char**);
