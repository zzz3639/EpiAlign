#Usage: julia this.jl index1 folder1 index2 folder2 stepsize p q outindex outfolder

#compress the long state sequence
function seq2sseq(seq)
    L = length(seq);
    Ans = Array(Any,2);
    s = seq[1];
    k = 1;
    t = 1;
    for i=2:L
        if seq[i]==s
	    k = k+1;
	else
	    s = seq[i];
	    t = t+1;
	    k = 1;
	end
    end
    Ans[1] = Array(Any,t);
    Ans[2] = Array(Any,t);
    s = seq[1];
    k = 1;
    t = 1;
    Ans[1][t] = s;
    for i=2:L
        if seq[i]==s
	    k = k+1;
	else
	    s = seq[i];
	    Ans[2][t] = k;
	    Ans[1][t+1] = s;
	    k = 1;
	    t = t+1;
	end
    end
    Ans[2][t] = k;
end


#read database1
f1 = open(ARGS[1],"r");
egindex1 = readlines(f1);
close(f1);
if egindex1[end]=="\n"
    L1 = length(egindex1) - 1;
else
    L1 = length(egindex1);
end

dbcontent1 = Array(Any,L1);

for i=1:L1
    if egindex1[i][end]=='\n'
        egindex1[i] = egindex1[i][1:end-1];
    end
    dbcontent1[i] = readdlm("$(ARGS[2])$(egindex1[i])",Int64);
end

#read database2
f2 = open(ARGS[3],"r");
egindex2 = readlines(f2);
close(f2);
if egindex2[end]=="\n"
    L2 = length(egindex2) - 1;
else
    L2 = length(egindex2);
end

dbcontent2 = Array(Any,L2);

for i=1:L2
    if egindex2[i][end]=='\n'
        egindex2[i] = egindex2[i][1:end-1];
    end
    dbcontent2[i] = readdlm("$(ARGS[4])$(egindex2[i])");
end

#initiate
chrnum = L1;
chrsize = Array(Int64,chrnum);
for i=1:chrnum
    chrsize[i] = 0;
end
for i=1:chrnum
    m = size(dbcontent1[i],1);
    for j=1:m
        chrsize[i] = chrsize[i] + dbcontent1[i][j,2];
    end
end
stepsize = parse(Float64,ARGS[5]);
p = parse(Float64,ARGS[6]);
q = parse(Float64,ARGS[7]);
thinit = q/(p+q);

#output index names
egindexmix = Array(Any,chrnum);
for i=1:chrnum
    egindexmix[i] = "$(egindex1[i]).mix";
end
f3 = open(ARGS[8],"w");
for i=1:chrnum
    write(f3,"$(egindexmix[i])\n");
end
close(f3);

#generate mixed epigenome chromosome by chromosome
srand();
for i=1:chrnum
    f4 = open("$(ARGS[9])$(egindexmix[i])","w");
    f5 = open("$(ARGS[9])$(egindexmix[i]).assign","w");
    seqthis1 = Array(Int64,chrsize[i]);
    seqthis2 = Array(Int64,chrsize[i]);
    k=0;
    for j=1:size(dbcontent1[i],1)
        seqthis1[k+1:k+dbcontent1[i][j,2]] = dbcontent1[i][j,1];
        k = k+dbcontent1[i][j,2];
    end
    k=0;
    for j=1:size(dbcontent2[i],1)
        seqthis2[k+1:k+dbcontent2[i][j,2]] = dbcontent2[i][j,1];
        k = k+dbcontent2[i][j,2];
    end
    seqthismix = Array(Int64,chrsize[i]);
    assignfrom = Array(Int64,chrsize[i]);
    v = Int64(ceil(1.0 * chrsize[i] / stepsize));
    k=0;
    l = k+1;
    u = min(k+stepsize,chrsize[i]);
    if rand()<thinit
        s = 1;
        assignfrom[l:u] = 1;
        seqthismix[l:u] = seqthis1[l:u];
    else
        s = 2;
        assignfrom[l:u] = 2;
	seqthismix[l:u] = seqthis2[l:u];
    end
    k = u;
    for j=2:v
        l = k+1;
	u = min(k+stepsize,chrsize[i]);
        if s==1
	    if rand()<p
	        s = 2;
		assignfrom[l:u] = 2;
		seqthismix[l:u] = seqthis2[l:u];
	    else
	        s = 1;
		assignfrom[l:u] = 1;
		seqthismix[l:u] = seqthis1[l:u];
	    end
	else
	    if rand()<q
	        s = 1;
		assignfrom[l:u] = 1;
		seqthismix[l:u] = seqthis1[l:u];
	    else
	        s = 2;
		assignfrom[l:u] = 2;
		seqthismix[l:u] = seqthis2[l:u];
	    end
	end
	k = u;
    end
    assignfromsseq = seq2sseq(assignfrom);
    seqthismixsseq = seq2sseq(seqthismix);

    Lassign = length(assignfromsseq[1]);
    Lmix = length(seqthismixsseq[1]);

    for j=1:Lassign
        write(f4,assignfromsseq[1][j]);
	write(f4," ");
	write(f4,assignfromsseq[2][j]);
	write(f4,"\n");
    end

    for j=1:Lmix
        write(f5,seqthismixsseq[1][j]);
	write(f5," ");
	write(f5,seqthismixsseq[2][j]);
	write(f5,"\n");
    end

    close(f4);
    close(f5);
end




