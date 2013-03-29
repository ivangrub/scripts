%Promoter Enrichment

clc

spacing = 250;

win = [-500 500];

idx = win/spacing;

Coverage = zeros(1,length(knownGene));


for j = 1:length(headers)
    x = eval(headers{j});
    fid = fopen(sprintf('HY_TSS_%s_Reads.txt',headers{j}),'w');
    fprintf(fid,'UCSC ID\tGene Abbr.\tGene Descript.\tChr\tStrand\tTSS\tReads\n');
    for i = 1:length(knownGene)

        if cell2mat(knownGene(i,5)) == 1
            TSS = round(cell2mat(knownGene(i,6))/spacing);
            TES = round(cell2mat(knownGene(i,7))/spacing);
            if TSS+idx(1) < 0
                Coverage(1,i) = sum(x.win.(knownGene{i,4})(2,1:TSS+idx(2)));
            else
                Coverage(1,i) = sum(x.win.(knownGene{i,4})(2,TSS+idx(1):TSS+idx(2)));
            end
        else
            TSS = round(cell2mat(knownGene(i,7))/spacing);
            TES = round(cell2mat(knownGene(i,6))/spacing);
            if TES-idx(2) < 0
                Coverage(1,i) = sum(x.win.(knownGene{i,4})(2,1:TSS));
            else
                Coverage(1,i) = sum(x.win.(knownGene{i,4})(2,TSS-idx(2):TSS-idx(1)));
            end
        end
        if cell2mat(knownGene(i,5)) == 1
            TSS = cell2mat(knownGene(i,6));
        else
            TSS = cell2mat(knownGene(i,7));
        end
        fprintf(fid,sprintf('%s\t%s\t%s\t%s\t%d\t%d\t%d\n',knownGene{i,1},knownGene{i,2}, ...
            knownGene{i,3},knownGene{i,4},cell2mat(knownGene(i,5)),TSS,Coverage(1,i)));
    end
    assignin('base',sprintf('%s_Coverage',headers{j}),Coverage);
    fclose(fid);clear fid
    
end

