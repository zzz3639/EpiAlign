function [] = MiddlePointSeg(pathmatch,chrom,pathseg,pathto,binsize)
%MiddlePointSeg('/home/zzz/EpiBLAST_baseline/test2/index', ...
%         '/home/zzz/StateDatabase/core_15/E003Files', ...
%         '/home/zzz/EpiBLAST_attention/output_attention/core_15/E003output/Segments/', ...
%         '/home/zzz/EpiBLAST_baseline/test2/',200)

chromnames = importdata(chrom);
L = length(chromnames);

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

mapchrom = containers.Map(chromstr',[1:L]);

D = importdata(pathmatch);
N = size(D,1);
strd = cell(2,1);
fout = fopen([pathto,'middle'],'w');
flength = fopen([pathto,'length'],'w');
for i=1:N
    strfull = D{i};
    k = findstr(strfull,' ');
    strd{1} = strfull(1:k-1);
    strd{2} = strfull(k+1:end);
    msseq = zeros(2,1);
    posa = zeros(2,1);
    posb = zeros(2,1);
    for j=1:2
        strthis = strd{j};
        k = findstr(strthis,'chr');
        if strthis(k+4)=='_'
            chromthis = strthis(k+3:k+3);
        else
            chromthis = strthis(k+3:k+4);
        end
        ka = findstr(strthis,'a');
        ka = ka(1);
        kb = findstr(strthis,'b');
        kc = findstr(strthis,'.anno');
        chrthis = mapchrom(chromthis);
        posa(j) = str2num(strthis(ka+1:kb-1));
        posb(j) = str2num(strthis(kb+1:kc-1));
        pathsseq = [pathseg,chromnames{chrthis},'/a',num2str(posa(j)),'b',num2str(posb(j)),'.sseq'];
        sseqthis = importdata(pathsseq);
        lsseq = size(sseqthis,1);
        if mod(lsseq,2)==0
            msseq(j) = sum(sseqthis(1:lsseq/2,2));
        else
            t = (lsseq-1)/2;
            msseq(j) = sum(sseqthis(1:t,2)) + floor(sseqthis(t+1,2)/2);
        end
        fprintf(flength,[num2str(size(sseqthis,1)),' ']);
    end
    fprintf(flength,'\n');
    fprintf(fout,[num2str(posa(1)*binsize + msseq(1)*binsize + 0.5*binsize),' ',num2str(posa(2)*binsize + msseq(2)*binsize + 0.5*binsize),'\n']);
end
fclose(fout);
fclose(flength);


end

