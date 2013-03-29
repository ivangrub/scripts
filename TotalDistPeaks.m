clc

headers = {'Oct4'};
comment = 'Lin';
load knownGene_mm9.mat
chr = importdata('mm9.chr.length');

for i = 1:length(headers)
    x = importdata(sprintf('%s_%s_Conden.txt',headers{i},comment));
    
    c = length(find(strcmp('chr1',knownGene(:,4))));
    TotalDistTSS = zeros(length(x.data),c);
    TotalDistTES = zeros(length(x.data),c);
    
    l = 1;
    for j= 1:length(chr.textdata)
        
        ch = find(strcmp(chr.textdata{j},knownGene(:,4)));
        p_ch = find(strcmp(x.textdata(2:end,5),chr.textdata{j}));
        
        for k = 1:length(ch)
            if cell2mat(knownGene(ch(k),5)) == 1
                TotalDistTSS(1:length(p_ch),l) = x.data(p_ch,1) - cell2mat(knownGene(ch(k),6));
                TotalDistTES(1:length(p_ch),l) = x.data(p_ch,1) - cell2mat(knownGene(ch(k),7));
            else
                TotalDistTSS(1:length(p_ch),l) = x.data(p_ch,1) - cell2mat(knownGene(ch(k),7));
                TotalDistTES(1:length(p_ch),l) = x.data(p_ch,1) - cell2mat(knownGene(ch(k),6));
            end
            l = l + 1; 
        end
    end
    
    
    assignin('base',sprintf('%s_Dist2TSSofPeaks',headers{i}),TotalDistTSS);
    assignin('base',sprintf('%s_Dist2TESofPeaks',headers{i}),TotalDistTES);
end