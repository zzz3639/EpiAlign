function [ SeqPerm, Fre ] = Perm1gram( Seq )
%PERM1GRAM Summary of this function goes here
%   Detailed explanation goes here
N=length(Seq);
Dic=unique(sort(Seq));
Fre=zeros(length(Dic),length(Dic));
for i=1:N-1
    c1=Seq(i);
    c2=Seq(i+1);
    v1=find(Dic==c1);
    v2=find(Dic==c2);
    Fre(v1,v2)=Fre(v1,v2)+1;
end
Fre=Fre./repmat(sum(Fre,2),1,length(Dic));

SeqPerm=Seq;
SeqPerm(1)=Seq(randi(N));

for i=2:N
    c=SeqPerm(i-1);
    v=find(Dic==c);
    vnew=mnrnd(1,Fre(v,:));
    SeqPerm(i)=Dic(find(vnew));
end

end

