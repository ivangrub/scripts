%Make the average reads plot for the TSS

load ESC_RNA_TSS.mm9.mat
load knownGene_mm9.mat
RNALimit = .01;

if isnan(RNALimit)
    knownGene = knownGene(isnan(GeneRNA(2,:)),:);
elseif RNALimit > 0.5
    knownGene = knownGene(GeneRNA(2,:) >= RNALimit,:);
else
    knownGene = knownGene(GeneRNA(2,:) <= RNALimit,:);
end

win = [-9000 -7000];
avgread = zeros((win+read_length-1)*2+1,length(headers)+1);
for i = 1:length(headers)
    x = eval(headers{i});
    gene = bp2TSS_win(x,knownGene,read_length,win);
end

