Type = {'Active' 'Inactive' 'NaN'};
RNALimit = [.01 .99 NaN];

load ESC_RNA_TSS.mm9.mat
load knownGene_mm9.mat

for i = 1:length(Type)

    fid = fopen(sprintf('%s_ESC_RNASeq.txt',Type{i}),'w');
   
    if isnan(RNALimit(i))
        known = knownGene(isnan(GeneRNA(2,1:49409)),:);
    elseif RNALimit(i) > 0.5
        known = knownGene(GeneRNA(2,1:49409) >= RNALimit(i),:);
    else
        known = knownGene(GeneRNA(2,1:49409) <= RNALimit(i),:);
    end
    
    for j = 1:length(known)
        fprintf(fid,sprintf('%s\t%s\t%s\t%s\t%d\t%d\t%d\n',char(known{j,2}),...
            char(known{j,2}),char(known{j,3}),char(known{j,4}),...
            cell2mat(known(j,5)),cell2mat(known(j,6)),cell2mat(known(j,7))));
    end
    fclose(fid);clear fid
end
