Reads = zeros(5500000,6);
A = zeros(5500000,2);
data = open('Francisco/ESC_TBP.mm9.mat');
data2 = open('Francisco/ESC_PolII_Ser5.mm9.mat');
data3 = open('Francisco/ESC_PolII_N.mm9.mat');
data4 = open('Francisco/ESC_TAF1.mm9.mat');
data5 = open('Francisco/ESC_PI.mm9.mat');
data6 = open('Francisco/ESC_MCPI.mm9.mat');
yy = strrep(fieldnames(data.chip),'chr','');
chr = cellstr(yy);
x = length(data.chip.(sprintf('chr%s',chr{1})));
tss1 = 1:20:x-40;
tss2 = 41:20:x;
y = length(tss1);
A(1:y,1:2) = [tss1' tss2'];
for j = 1:y
    Reads(j,1) = sum(data.chip.chr1(A(j,1):A(j,2)));
    Reads(j,2) = sum(data2.chip.chr1(A(j,1):A(j,2)));
    Reads(j,3) = sum(data3.chip.chr1(A(j,1):A(j,2)));
    Reads(j,4) = sum(data4.chip.chr1(A(j,1):A(j,2)));
    Reads(j,5) = sum(data5.chip.chr1(A(j,1):A(j,2)));
    Reads(j,6) = sum(data6.chip.chr1(A(j,1):A(j,2)));
end
for i = 2:length(chr);
    x = length(data.chip.(sprintf('chr%s',chr{i})));
    tss1 = 1:20:x-40;
    tss2 = 41:20:x;
    a = length(tss1);
    C = [tss1' tss2'];
    A(y:y+a,1:2) = C;
    for j = 1:y
       NewReads(:,1) = sum(data.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
       NewReads(:,2) = sum(data2.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
       NewReads(:,3) = sum(data3.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
       NewReads(:,4) = sum(data4.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
       NewReads(:,5) = sum(data5.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
       NewReads(:,6) = sum(data6.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
       Reads = vertcat(Reads,NewReads);
    end
end

FullGenome_Windows = A;

chr_num = yy;
%Note: This code does not save the chromosome number associated with each
%TSS window. Use chr_save.m to create a cellular vector with those values.

clear A C y tss1 tss2  data* yy 
