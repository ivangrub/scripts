function [chipout] ...
    = importchipseq(data,CpG,knownGene)

[m,n] = size(data);

%Open enrichment data matrix
chip_raw = data;
%chipFR_raw = data(:,[2:2:n]);

%Normalize with respect to the background
chip = normalize(chip_raw);
%chipFR = normalize(chipFR_raw);

%Sort the data into those with CpG islands and the remaining data for
%both ChIP and ChIP-FR
%[wcpggene wcpg wocpggene wocpg] = sortout(chip(:,1:2),CpG,knownGene);
%[uniquegeneFR uniquewcpgFR wcpggeneFR wcpgFR wocpggeneFR wocpgFR] = sortout(chipFR(:,1:3),CpG,knownGene);

chipout = chip(:,1:3);
%chipoutFR = chipFR(:,1:3);
end