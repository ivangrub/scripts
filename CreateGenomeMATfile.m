clear, clc

genomename = 'dm3';

data = importdata(sprintf('%s.chr.length',genomename));
len = data.data;
chr = data.textdata;

for i = 1:length(chr)
    genome.(chr{i}) = fastaread(sprintf('/Users/ivang/Documents/MATLAB/Tjian Lab/ChIP-Seq/DM3CHR/%s.fa',chr{i}));
end

save(sprintf('genome_%s.mat',genomename),'-v7.3')
