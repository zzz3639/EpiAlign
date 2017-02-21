#ifndef CustomFunction_H
#define CustomFunction_H

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"WatermanFun.h"
#include"StateIO.h"


void Custom_Init(unsigned char**,unsigned short**,float**,int*,unsigned char**,unsigned short**,float**,int*,float*,void**,const char*);


void Custom_Free(void*);


float Custom_MatchingFunction(unsigned char, unsigned char, int, int, void*);


float Custom_GapFunction(float, int, char, void*);


#endif
