%Make the average reads plot for the TSS

clc

win = [-9000 -7000];                        %Relative to TSS (i.e. negative sign means upstream)
IP = 1;

load ESC_RNA_TSS.mm9.mat
load knownGene_mm9.mat

RNALimit = .99;

if isnan(RNALimit)
    known = knownGene(isnan(GeneRNA(2,1:49409)),:);
elseif RNALimit > 0.5
    known = knownGene(GeneRNA(2,1:49409) >= RNALimit,:);
else
    known = knownGene(GeneRNA(2,1:49409) <= RNALimit,:);
end

x = eval(headers{IP});
datawin = x.win.chr1(1,2)-x.win.chr1(1,1);

gene = cell(1,length(known));
coverage = zeros(1,length(known));

for i = 1:length(known)
    gene(1,i) = known(i,2);
    if cell2mat(known(i,5)) == 1
        win2 = [cell2mat(known(i,6))+win(1) cell2mat(known(i,6))+win(2)];
    else 
        win2 = [cell2mat(known(i,7))-win(1) cell2mat(known(i,6))-win(2)];
    end
    if win2(1) > 0 && win2(2) > 0
        coverage(1,i) = sum(x.win.(known{i,4})(2,floor((win2(1))/datawin):floor((win2(2)-1)/datawin)));
    end
end

[~,j] = sort(coverage,'descend');
sortedgene = gene(j);
sorted = coverage(j);

clear avg x i