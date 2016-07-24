function [] = putchrname( Fileto, Names )
%PUTCHRNAME Summary of this function goes here
%   Detailed explanation goes here
fp = fopen(Fileto,'w');
N = length(Names);
for i=1:N
    fprintf(fp,'%s\n',Names{i});
end
fclose(fp);
end

