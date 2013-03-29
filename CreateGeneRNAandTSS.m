clc
load knownGene_mm9.mat
load FPKM_ESC.mat

TSS = cell(length(knownGene),5);

for i = 1:length(knownGene)
    if cell2mat(knownGene(i,5)) == 1
        TSS(i,:) = {knownGene{i,1} knownGene{i,4} cell2mat(knownGene(i,5)) cell2mat(knownGene(i,6)) cell2mat(knownGene(i,7))};
    else
        TSS(i,:) = {knownGene{i,1} knownGene{i,4} cell2mat(knownGene(i,5)) cell2mat(knownGene(i,7)) cell2mat(knownGene(i,6))};
    end
end

GeneRNA = zeros(2,length(TSS));
for k = 1:length(TSS)
    LL = find(strcmp(TSS(k,1),FPKM_ESC(:,1)));
    if ~isempty(LL) 
        GeneRNA(1,k) = max(FPKM_ESC{LL,4});
        GeneRNA(2,k) = min(FPKM_ESC{LL,5});
    else 
        GeneRNA(1:2,k) = [0;1];
    end
end
