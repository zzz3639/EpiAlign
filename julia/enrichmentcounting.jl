#Usage: julia this.jl genelist regionlist
#ARGS = Array(UTF8String,2);
#ARGS[1] = "/home/zzz/EpiEnrich/escellgene"
#ARGS[2] = "/home/zzz/EpiBLAST_baseline/test2/position"
#ARGS[2] = "/home/zzz/EpiBLAST_attention/output_attention/core_15/E003output/test/position"

#read genelist
fingene = open(ARGS[1],"r");
genelist = readlines(fingene);
close(fingene);
if genelist[end]=="\n"
    L1 = length(genelist) - 1;
else
    L1 = length(genelist);
end
genearray = Array(Any,L1,7);

for i=1:L1
    genearray[i,:] = split(genelist[i][1:end-1],'\t');
end

for i=1:L1
    genearray[i,6] = parse(Int64,genearray[i,6]);
    genearray[i,7] = parse(Int64,genearray[i,7]);
end

#read regionlist
finregion = open(ARGS[2],"r");
regionlist = readlines(finregion);
close(finregion);
if regionlist[end]=="\n"
    L2 = length(regionlist) - 1;
else
    L2 = length(regionlist);
end
regionarray = Array(Any,L2,3);

for i=1:L2
    strsplit = split(regionlist[i],' ');
    parsestr = split(strsplit[1],':');
    regionarray[i,1] = parsestr[1];
    parsestr = split(parsestr[2],'-');
    regionarray[i,2] = parse(Int64,parsestr[1]);
    regionarray[i,3] = parse(Int64,parsestr[2]);
end

#count overlap number
nol = 0;
for i=1:L1
    for j=1:L2
        if regionarray[j,1]!=genearray[i,4]
            continue;
	end
	if regionarray[j,3] < genearray[i,6]
	    continue;
	end
	if regionarray[j,2] > genearray[i,7]
	    continue;
	end
	nol += 1;
	println("gene:$(genearray[i,2]):$(genearray[i,4]):$(genearray[i,6])-$(genearray[i,7]) region:$(regionarray[j,1]):$(regionarray[j,2])-$(regionarray[j,3])\n");
    end
end

println("\n$nol\n")

