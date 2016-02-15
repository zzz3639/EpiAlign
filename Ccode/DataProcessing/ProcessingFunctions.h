#ifndef PROCESSINGFUNCTIONS_H
#define PROCESSINGFUNCTIONS_H

#include<stdio.h>
#include<string>
#include<list>
#include<utility>
namespace CONVERT
{
    void Seq2SSeq(std::list<short>, std::list<pair<short,int> >&);
}

namespace READWRITE
{
    void ReadSequence(string, std::list<short>&);
}


#endif