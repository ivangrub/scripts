%Sorts the raw enrichment data into different subpopulations that can then
%be used in the PCA

function [uniquegene uniquewcpg wocpggene wocpg] = sortout(rawdata,CpG,knownGene)

[m,n] = size(rawdata);
j = 1;
k = 1;
l = 1;
for i = 1:length(rawdata)
    if CpG(i,1) ~=0 && knownGene{i,10} == 1
        uniquewcpg(j,1:n) = rawdata(i,1:n);
        uniquegene{j,1} = knownGene{i,2};
        uniquegene{j,2} = knownGene{i,3};
        j = j + 1;
%     elseif CpG(i,1) ~= 0 && knownGene{i,10} == 0
%         wcpg(k,1:n) = rawdata(i,1:n);
%         wcpggene{k,1} = knownGene{i,2};
%         wcpggene{k,2} = knownGene{i,3};
%         k = k + 1;
    else
        wocpg(l,1:n) = rawdata(i,1:n);
        wocpggene{l,1} = knownGene{i,2};
        wocpggene{l,2} = knownGene{i,3};
        l = l + 1;
    end
end
end