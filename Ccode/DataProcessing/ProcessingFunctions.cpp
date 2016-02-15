#include"ProcessingFunctions.h"

namespace CONVERT
{
	void Seq2SSeq(std::list<short> InSeq, list<pair<short,int> > &OutSeq)
	{
		std:list<short>::iterator itin;
		std:list<pair<short,int> >::iterator itout;
		/*Empty OutSeq*/
		while(!itout.empty()){
            itout.pop_front();
		}
		/*Assign OutSeq*/
		short state=-1;
		std::pair<short,int> temppair(-1,0);
		for(itin=InSeq.begin();itin!=InSeq.end();itin++){
            if(state!=*it){
            	if(temppair.first!=-1)
            		OutSeq.push_back(temppair);
                temppair=std::pair<short,int>(*it,1);
                state=*it;
            }
            else{
            	temppair.second++;
            }
		}
		if(temp.first!=-1)
			OutSeq.push_back(temppair)
	}
}

namespace READWRITE
{
	void ReadSequence(string, std::list<short>&)
	{
		return;
	}
}