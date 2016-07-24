function [StateSeq, L] = EpiBrowser(Pos,Opt)
n = Opt.n;
chr = Pos(1);
pos = floor((Pos(2)-1)/Opt.w);
a = max(floor(pos-Opt.m/2),0);
a = floor(a/Opt.k)*Opt.k;
b = a+Opt.m;
L = pos-a;
FileIndex = cell(n,1);
StateSeq = cell(n,1);
for i=1:n
    FileIndex{i} = importdata(Opt.index{i});
end

for i=1:n
    Seqthis = importdata([Opt.folder{i},'/',FileIndex{i}{chr},'/a',num2str(a),'b',num2str(b),'.seq']);
    l = length(Seqthis);
    strthis = '';
    for j=1:l
        strthis = [strthis,Seqthis{j}];
    end
    StateSeq{i} = strthis;
end
    
end