function [ Top, Map ] = CountComb( StateSeq, l, topn )
%COUNTFRE Summary of this function goes here
%   Detailed explanation goes here
if size(StateSeq,1)~=1
    StateSeq=StateSeq';
end
N=length(StateSeq);
StateSeq=char(StateSeq+'0');
MapKey=cell(1,1);
MapKey{1}=StateSeq(1:l);
MapValue=0;
Map=containers.Map(MapKey,MapValue);

for i=1:N-l+1
    KeyThis=StateSeq(i:i+l-1);
    if isKey(Map,KeyThis)
        Map(KeyThis)=Map(KeyThis)+1;
    else
        Map(KeyThis)=1;
    end
end

TopCount=zeros(topn,1);
TopComb=zeros(topn,l);

for KeyThis=Map.keys
    if Map(KeyThis{1})>TopCount(1)
        TopCount(1)=Map(KeyThis{1});
        TopComb(1,:)=KeyThis{1}-'0';
    end
    [u,v]=sort(TopCount);
    TopCount=TopCount(v,1);
    TopComb=TopComb(v,:);
end

Top.TopCount=TopCount;
Top.TopComb=TopComb;

end

