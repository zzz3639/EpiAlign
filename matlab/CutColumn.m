function [] = CutColumn( pathf, patht, charsp )
%CUTCOLUMN Summary of this function goes here
%   Detailed explanation goes here
Lines = importdata(pathf);
L = length(Lines);

fout = fopen(patht,'w');
for i=1:L
    strthis = Lines{i};
    k = strfind(strthis,charsp);
    fprintf(fout,[strthis(1:k-1),'\n']);
end
fclose(fout);

end

