%Make the average reads plot for the TSS

load ESC_RNA_TSS.mm9.mat
load knownGene.mm9.mat
RNALimit = NaN;

if isnan(RNALimit)
    knownGene = knownGene(isnan(GeneRNA(2,:)),:);
elseif RNALimit > 0.5
    knownGene = knownGene(GeneRNA(2,:) >= RNALimit,:);
else
    knownGene = knownGene(GeneRNA(2,:) <= RNALimit,:);
end

win = [-1000 -900];
avgread = zeros((win+read_length-1)*2+1,length(headers)+1);
for i = 1:length(headers)
    x = eval(headers{i});
    avg = bp2TSS_wRNA(x,knownGene,read_length,win);
    avgread(:,[1 i+1]) = avg;
end

clear avg x i