function [] = CopyMaxAnnoBase( chrom, pathf, patht, step, win, binsize, numrecord)
% Usage: CopyMaxAnnoBase( '/home/zzz/StateDatabase/core_15/E003Files', ...
%        '/home/zzz/EpiBLAST_baseline/test/', '/home/zzz/EpiBLAST_baseline/test2/', ...
%        500, 2500, 200, 100)
%   Detailed explanation goes here
chromnames = importdata(chrom);
L = length(chromnames);
Mfull = zeros(0,23);
for i=1:L
    M = importdata([pathf,chromnames{i},'.baselinescore']);
    if size(M,1)==0
        continue;
    end
    M = [M,i*ones(size(M,1),1),[0:size(M,1)-1]'];
    Mfull = [Mfull;M];
end

A = Mfull(:,3) - Mfull(:,1);
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
    system(['cp ',pathf,chromnames{Mtop(end-i,22)},'/a',num2str(Mtop(end-i,23)*step),'b',num2str(Mtop(end-i,23)*step+win),'.sseq.anno ',patht,'chr',chromstr{Mtop(end-i,22)},'_a',num2str(Mtop(end-i,23)*step),'b',num2str(Mtop(end-i,23)*step+win),'.anno']);
end

fidx = fopen([patht,'index'],'w');
fpos = fopen([patht,'position'],'w');
flist = fopen([patht,'tablelist'],'w');
for i=0:numrecord-1
    annolines = importdata([patht,'chr',chromstr{Mtop(end-i,22)},'_a',num2str(Mtop(end-i,23)*step),'b',num2str(Mtop(end-i,23)*step+win),'.anno']);
    matchregion = annolines.textdata{3,2};
    k = findstr(matchregion,'chr');
    if matchregion(k+4)==' '
        chromthis = matchregion(k+3:k+3);
    else
        chromthis = matchregion(k+3:k+4);
    end
    ka = findstr(matchregion,'a');
    kb = findstr(matchregion,'b');
    posthis = str2num(matchregion(ka+1:kb-1));
    if posthis==Mtop(end-i,23)*step&&strcmp(chromthis,chromstr{Mtop(end-i,22)})
        matchregion = annolines.textdata{2,2};
        k = findstr(matchregion,'chr');
        if matchregion(k+4)==' '
            chromthis = matchregion(k+3:k+3);
        else
            chromthis = matchregion(k+3:k+4);
        end
        ka = findstr(matchregion,'a');
        kb = findstr(matchregion,'b');
        posthis = str2num(matchregion(ka+1:kb-1));
    end
    fprintf(fpos,['chr',chromstr{Mtop(end-i,22)},':',num2str(Mtop(end-i,23)*step*binsize),'-',num2str((Mtop(end-i,23)*step+win)*binsize),' ']);
    fprintf(fpos,['chr',chromthis,':',num2str(posthis*binsize),'-',num2str((posthis+win)*binsize),'\n']);
    fprintf(fidx,['chr',chromstr{Mtop(end-i,22)},'_a',num2str(Mtop(end-i,23)*step),'b',num2str(Mtop(end-i,23)*step+win),'.anno ']);
    fprintf(fidx,['chr',chromthis,'_a',num2str(posthis),'b',num2str(posthis+win),'.anno\n']);
    fprintf(flist,['chr',chromstr{Mtop(end-i,22)},':',num2str(Mtop(end-i,23)*step*binsize),'-',num2str((Mtop(end-i,23)*step+win)*binsize),'\n']);
    fprintf(flist,['chr',chromthis,':',num2str(posthis*binsize),'-',num2str((posthis+win)*binsize),'\n']);
end
fclose(fpos);
fclose(fidx);
fclose(flist);


end

