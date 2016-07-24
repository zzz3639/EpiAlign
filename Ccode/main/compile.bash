gcc TwoRegions.c WatermanFun.c CustomFunction.c StateIO.c -o ../bin/TwoRegions.out
gcc cut_sseq.c WatermanFun.c StateIO.c -o ../bin/cut_sseq.out
gcc DatabaseSearch.c WatermanFun.c StateIO.c CustomFunction.c -g -o ../bin/DatabaseSearch.out
