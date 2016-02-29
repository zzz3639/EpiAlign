function [] = stat( Data, Outfilename, MotifLength )
%STAT Summary of this function goes here
%   Detailed explanation goes here

%const
Ntop=200;
Zoom=20;
% initialize
Seq=Data(:,1);
Len=Data(:,2);
Dic=unique(sort(Seq));
Quies=Dic(end);

%count segment lengths
Ave=mean(Len);
v=find(Seq~=Quies);
AveNoQuies=mean(Len(v));

%count frequencies of short combinations
[ Top, Map ] = CountComb( Seq, MotifLength, Ntop );
[ TopPerm, MapPerm ] = CountComb( Rperm(Seq), MotifLength, Ntop );

%count combination baseline
[Gram1Seq, TransMatrix]=Perm1gram(Seq);
imshow(kron(TransMatrix,ones(Zoom,Zoom)));
[ TopGram1, MapGram1 ] = CountComb( Gram1Seq, MotifLength, Ntop );

FoldMotifs=zeros(Ntop,1);
for i=1:Ntop
    strthis=char(Top.TopComb(i,:)+'0');
    countthis=Top.TopCount(i)+1;
    if isKey(MapGram1,strthis)
        countgram1=MapGram1(strthis)+1;
    else
        countgram1=1;
    end
    FoldMotifs(i)=countthis/countgram1;
end

save([Outfilename,'.mat']);
WriteHTMLFile(Seq,Len,Outfilename);

end

function []=WriteHTMLFile(Seq,Len,Outfilename)
N=length(Seq);
SeqStr=char(Seq'+'0');
fid = fopen([Outfilename,'.html'],'wt');
fprintf(fid, '<html>\n <body>\n');
for i=1:N
    fprintf(fid,['<abbr title="',num2str(Len(i)),'">',SeqStr(i),'</abbr>\n']);
end
fprintf(fid, '</body> </html>\n');
fclose(fid);
end



