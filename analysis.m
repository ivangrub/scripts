%Organize the analysis of the ChIP-Seq data

%Run the principle component analysis and normalize the data to the
%background.

ChipSeqOptions;

%Spread out the data along the log axis

LoggedCell = log2(1+CellType);

%Imports the sorted data that only identifies regions that have positive
%enrichments.
if Log == 1
    [cpgchip] = importchipseq(LoggedCell,CpG,knownGene);
else
    [cpgchip] = importchipseq(CellType,CpG,knownGene);
end

[k,l] = size(cpgchip);

%Run the PCA for each one individually
s1 = 'cpgchip';
s2 = 'restchip';
i = 1;
for n = 1
    x = eval(sprintf('s%d',n));
    [PC Zscore coeff latent cumvar] = PCA(eval(x));
    Cumvar(n,:) = cumvar;
    EigTotal(:,n) = latent;
    i = i + l;
    assignin('base',sprintf('%sPCTotal%d',option,n),PC(:,1:l));
    assignin('base',sprintf('%sZscore%d',option,n),Zscore(:,1:l));
    assignin('base',sprintf('%sCoeff%d',option,n),coeff);
end
clear s*


