function [ filenames_sort ] = sortchrname( filenames )
%SORTNAME Summary of this function goes here
%   Detailed explanation goes here
N = length(filenames);
chrnum = zeros(N,1);
for i=1:N
    str = filenames{i};
    k = findstr(str,'chr');
    if str(k+3)>='0'&&str(k+3)<='9'
        if str(k+4)=='_'
            chrnum(i) = str2num(str(k+3:k+3));
        else
            chrnum(i) = str2num(str(k+3:k+4));
        end
    else
        chrnum(i) = str(k+3) + N;
    end
end

[u,v] = sort(chrnum);
filenames_sort = filenames(v,:);

end

