%Run this script when you approximately know the length of the full genome
%analysis. New_coverage_wholegenome.m will be easier for first time use
%since it just conacates the arrays, but it does take longer.

clear

headers = {'PreTAF7L' 'PrePI' 'PostTAF7L' 'PostPI'};
data = open('Haiying/PreDiff/HY_TAF7L.mm9.mat');
data2 = open('Haiying/PreDiff/HY_PI.mm9.mat');
data3 = open('Haiying/PostDiff/HY_TAF7L.mm9.mat');
data4 = open('Haiying/PostDiff/HY_PI.mm9.mat');
chr_save;
m = length(Chr);
Reads = zeros(m,4);
A = zeros(m,2);

yy = strrep(fieldnames(data.chip),'chr','');
chr = cellstr(yy);
y = 1;
k = 1;
for i = 1:length(chr);
    x = length(data.chip.(sprintf('chr%s',chr{i})));
    tss1 = 1:5:x-10;
    tss2 = 11:5:x;
    a = length(tss1);
    C = [tss1' tss2'];
    A(y:y+a-1,1:2) = C;
    for j = 1:a
       NewReads(1,1) = sum(data.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
       NewReads(1,2) = sum(data2.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
       NewReads(1,3) = sum(data3.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
       NewReads(1,4) = sum(data4.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
       Reads(k,1:4) = NewReads;
       k = k + 1;
    end
    y = y + a;
end

RealLength = y - 1;
FullGenome_Windows = A(1:RealLength,:);
Reads = Reads(1:RealLength,:);
chr_num = yy;

%Note: This code does not save the chromosome number associated with each
%TSS window. Use chr_save.m to create a cellular vector with those values.

clear A C y tss1 tss2  data* yy NewReads a k chr i j x m


