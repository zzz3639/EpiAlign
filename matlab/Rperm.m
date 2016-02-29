function [ SeqPerm ] = Rperm( Seq )
%RPERM Summary of this function goes here
%   Detailed explanation goes here
N=length(Seq);
SeqPerm=Seq;
SeqPerm(1)=Seq(randi(N));
Dic=unique(sort(Seq));
Fre=cell(length(Dic),1);
for i=1:length(Dic)
    Fre{i}=find(Seq~=Dic(i));
end
for i=2:N
    i;
    c=SeqPerm(i-1);
    v=Fre{find(Dic==c)};
    SeqPerm(i)=Seq(v(randi(length(v))));
end

end

