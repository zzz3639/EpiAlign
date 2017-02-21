function [] = SegCopy( pathf,patht,v,a,l )
%SEGCOPY Summary of this function goes here
%   Detailed explanation goes here
n = length(v);
for i=1:n
    k = v(i);
    system(['cp ',pathf,'a',int2str((k-1)*l),'b',int2str((k-1)*l+a),'* ',patht]);
end

end

