gcc TwoRegions.c WatermanFun.c CustomFunction.c StateIO.c -O3 -o ../bin/TwoRegions.out
gcc cut_sseq.c WatermanFun.c StateIO.c -O3 -o ../bin/cut_sseq.out
gcc DatabaseSearch.c WatermanFun.c StateIO.c CustomFunction.c -O3 -o ../bin/DatabaseSearch.out
gcc DatabaseSearchNative.c WatermanFun.c StateIO.c CustomFunction.c -O3 -o ../bin/DatabaseSearchNative.out
gcc Algn2Ano.c StateIO.c -O3 -o ../bin/Algn2Ano.out
g++-5 -std=c++11 FakeChromosomeGenerator.cpp WatermanFun.c StateIO.c -O3 -o ../bin/FakeChromosomeGenerator.out
gcc -O3 DatabaseSearch_baseline1.c WatermanFun.c StateIO.c -lm -o ../bin/DatabaseSearch_baseline.out
