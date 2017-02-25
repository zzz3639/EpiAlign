function [] = CopyMaxAnno(chrom, pathf, patht, step, win, binsize, numrecord)
% CopyMaxAnno('/home/zzz/StateDatabase/core_15/E003Files', ...
%             '/home/zzz/EpiBLAST_attention/output_attention/core_15/E003output/', ...
%             '/home/zzz/EpiBLAST_attention/output_attention/core_15/E003output/test/', 500, 2500, 200, 100);
%   Detailed explanation goes here

chromnames = importdata(chrom);
L = length(chromnames);
Mfull = zeros(0,7);
for i=1:L
    M = importdata([pathf,chromnames{i},'.score']);
    if size(M,1)==0
        continue;
    end
    M = [M,i*ones(size(M,1),1),[0:size(M,1)-1]'];
    Mfull = [Mfull;M];
end

A = (Mfull(:,4) - Mfull(:,2))./Mfull(:,2);
[u,v] = sort(A);

Mtop = Mfull(v(end-numrecord+1:end),:);


chromstr = cell(L,1);
for i=1:L
    str = chromnames{i};
    k = findstr(str,'chr');
    if str(k+3)>='0'&&str(k+3)<='9'
        if str(k+4)=='_'
            chromstr{i} = str(k+3:k+3);
        else
            chromstr{i} = str(k+3:k+4);
        end
    else
        chromstr{i} = str(k+3);
    end
end

for i=0:numrecord-1
    system(['cp ',pathf,chromnames{Mtop(end-i,6)},'/a',num2str(Mtop(end-i,7)*step),'b',num2str(Mtop(end-i,7)*step+win),'.sseq.anno ',patht,'chr',chromstr{Mtop(end-i,6)},'_a',num2str(Mtop(end-i,7)*step),'b',num2str(Mtop(end-i,7)*step+win),'.anno']);
end

flength = fopen([patht,'length'],'w');
for i=0:numrecord-1
    fprintf(flength,[num2str(Mtop(end-i,1)),'\n']);
end
fclose(flength);

fidx = fopen([patht,'index'],'w');
fpos = fopen([patht,'position'],'w');
flist = fopen([patht,'tablelist'],'w');

for i=0:numrecord-1
    fprintf(fidx,['chr',chromstr{Mtop(end-i,6)},'_a',num2str(Mtop(end-i,7)*step),'b',num2str(Mtop(end-i,7)*step+win),'.anno\n']);
    C = importdata([patht,'chr',chromstr{Mtop(end-i,6)},'_a',num2str(Mtop(end-i,7)*step),'b',num2str(Mtop(end-i,7)*step+win),'.anno'],'\n');
    strthis = C{9};
    ks = strfind(strthis,' ');
    kq = strfind(strthis,'Query:');
    kd = strfind(strthis,'Database:');
    fprintf(fpos,strthis(kq+6:ks));
    fprintf(fpos,[strthis(kd+9:end),'\n']);
    fprintf(flist,[strthis(kq+6:ks),'\n']);
    fprintf(flist,[strthis(kd+9:end),'\n']);
end

fclose(fidx);
fclose(fpos);
fclose(flist);

end

