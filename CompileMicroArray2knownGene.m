%Microarray Gene and Data into knownGene table

load knownGene.mm9.mat
load ../../RNA-Seq/known2Ref.mm9.mat

CompiledMicroData = cell(length(knownGene),5);

for i = 1:length(knownGene)
    x = find(strcmp(knownGene(i,1),known2Ref(:,1)));
    if ~isempty(x)
        y = find(strcmp(known2Ref(x,2),MicroArrayGenes(:,1)));
        CompiledMicroData(i,:) = {knownGene(i,1) known2Ref(x,2) knownGene(i,2) ...
            MicroArrayData(2,x) MicroArrayData(1,x)};
    else
        CompiledMicroData(i,:) = {NaN NaN NaN NaN NaN};
    end
end