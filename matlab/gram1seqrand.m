function [seqgram1randfull, seqgram1randfull_sseq, seqgram1randfull_num] = gram1seqrand(testfile)
seq = testfile(:,1);
m = sum(testfile(:,2));
seqfull = zeros(m,1);
n=0;
for i=1:length(seq)
    seqfull(n+1:n+testfile(i,2),1) = seq(i);
    n = n+testfile(i,2);
end
seqgram1randfull = Perm1gram(seqfull,length(seqfull));

D = [seqgram1randfull;max(seqfull)+1]-[0;seqgram1randfull];
v = find(D~=0);
seqgram1randfull_num = v(2:end)-v(1:end-1);
seqgram1randfull_sseq = seqgram1randfull(v(1:end-1));

end