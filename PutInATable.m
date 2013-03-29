%Put the condensed/significant genes into a table

clc

%% Import the GeneID and set controls/comments

A = importdata('GeneID');
GeneID = A.data;
UCSC = A.textdata(2:end,1);

control = find(strcmp('PI',headers));
comment = 'CC_ESC_Lin';

%% Convert to GeneID and add chr into knownGene

GID = NaN(length(knownGene),1);
k = 1;
for i = 1:length(knownGene)
    a = find(strcmp(knownGene{i,1},UCSC));
    if isempty(a) == 0
        GID(k,1) = GeneID(a);
        k = k + 1;
    else
        k = k + 1;
    end
end
knownGene(:,4) = strcat('chr',knownGene(:,4));

%% Initialize the necessary matrices and options

Location = {'TSS' 'TSS' 'TES' 'TES'};
Orientation = {'Upstream' 'Downstream'};
Columns = [3 5 8 10];
Win = [1 1 6 6];
SigTable = cell(length(knownGene),17);
CondenTable = cell(length(knownGene),17);

%% Run through the significant and condensed structures

l = 1;
n = 1;
for i = 1:(control-1)
    
    sig = eval(sprintf('%s_significant',headers{i}));
    
    %Printing table for significant subset
    fprintf('Printing table for %s significant\n',headers{i})
    fid = fopen(sprintf('%s_%s_SigTable.txt',headers{i},comment),'w');
    fprintf(fid,'UCSC-ID\t GENE-ID\t Gene Abbreviation\t Gene Descript.\t Chr\t Strand\t TSS\t TES\t Feature\t Relative Position\t Inter or Intra\t Middle of Window\t Distance to Feature\t Left Edge\t Right Edge\t Reads\t Enrichment\n');
    for j = 1:length(knownGene)
        a = find((j == sig.Distance(:,2)));
        if isempty(a)
            continue
        end
        for k = 1:length(a)
            if knownGene{j,5} == 1 && strcmp(Location(b(k)),'TSS')
                if sign(sig.Distance(a(k),Columns(b(k)))) == -1
                    o = 1;
                else
                    o = 2;
                end
            elseif knownGene{j,5} == -1 && strcmp(Location(b(k)),'TSS')
                if sign(sig.Distance(a(k),Columns(b(k)))) == -1
                    o = 2;
                else
                    o = 1;
                end
            elseif knownGene{j,5} == 1 && strcmp(Location(b(k)),'TES')
                if sign(sig.Distance(a(k),Columns(b(k)))) == -1
                    o = 1;
                else
                    o = 2;
                end
            elseif knownGene{j,5} == -1 && strcmp(Location(b(k)), 'TES')
                if sign(sig.Distance(a(k),Columns(b(k)))) == -1
                    o = 2;
                else
                    o = 1;
                end
            end
            if (sig.Distance(a(k),Win(b(k))) > knownGene{j,6}) && ...
                    (sig.Distance(a(k),Win(b(k))) < knownGene{j,7})
                SigTable(l,1:17) = {knownGene(j,1) GID(j) knownGene{j,2:7} ...
                    Location(b(k)) Orientation(o) 'Intragenic' sig.Distance(a(k),Win(b(k))) knownGene{j,5}*sig.Distance(a(k),Columns(b(k))) ...
                    sig.Window(a(k),1) sig.Window(a(k),2) sig.Reads(a(k),1) sig.Enrichment(a(k),1)};
            else
                SigTable(l,1:17) = {knownGene(j,1) GID(j) knownGene{j,2:7} ...
                    Location(b(k)) Orientation(o) 'Intergenic' sig.Distance(a(k),Win(b(k))) knownGene{j,5}*sig.Distance(a(k),Columns(b(k))) ...
                    sig.Window(a(k),1) sig.Window(a(k),2) sig.Reads(a(k),1) sig.Enrichment(a(k),1)};
            end
            fprintf(fid,sprintf('%s\t %d\t %s\t %s\t %s\t %d\t %d\t %d\t %s\t %s\t %s\t %d\t %d\t %d\t %d\t %d\t %d\n', ...
                char(SigTable{l,1}),cell2mat(SigTable(l,2)),char(SigTable{l,3}),char(SigTable{l,4}), ...
                char(SigTable{l,5}),cell2mat(SigTable(l,6)),cell2mat(SigTable(l,7)), ...
                cell2mat(SigTable(l,8)),char(SigTable{l,9}),char(SigTable{l,10}), ...
                char(SigTable{l,11}),cell2mat(SigTable(l,12)),cell2mat(SigTable(l,13)), ...
                cell2mat(SigTable(l,14)),cell2mat(SigTable(l,15)),cell2mat(SigTable(l,16)),cell2mat(SigTable(l,17))));
            l = l + 1;
        end
    end
    if l < length(knownGene)
        SigTable = SigTable(1:l-1,:);
    end
    fclose(fid);clear fid
    
    fprintf('Printing table for %s condensed\n',headers{i})
    %Printing table for condensed windows
    fid = fopen(sprintf('%s_%s_CondenTable.txt',headers{i},comment),'w');
    fprintf(fid,'UCSC-ID\t GENE-ID\t Gene Abbreviation\t Gene Descript.\t Chr\t Strand\t TSS\t TES\t Feature\t Relative Position\t Inter or Intra\t Middle of Window\t Distance to Feature\t Left Edge\t Right Edge\t Reads\t Enrichment\n');
    for j = 1:length(knownGene)
            conden = eval(sprintf('%s_condensed',headers{i}));
            a = find((j == conden.Distance(:,2)));
            if isempty(a)
                continue
            end
            
            for m = 1:length(a)
                if knownGene{j,5} == 1 && (k == 1 || k == 2)
                    if sign(conden.Distance(a(m),3)) == -1
                        o = 1;
                    else
                        o = 2;
                    end
                elseif knownGene{j,5} == -1 && (k == 1 || k == 2)
                    if sign(conden.Distance(a(m),3)) == -1
                        o = 2;
                    else
                        o = 1;
                    end
                elseif knownGene{j,5} == 1 && (k == 3 || k == 4)
                    if sign(conden.Distance(a(m),3)) == -1
                        o = 1;
                    else
                        o = 2;
                    end
                elseif knownGene{j,5} == -1 && (k == 3 || k == 4)
                    if sign(conden.Distance(a(m),3)) == -1
                        o = 2;
                    else
                        o = 1;
                    end
                end
                if (conden.Distance(a(m),1) > knownGene{j,6}) && ...
                        (conden.Distance(a(m),1) < knownGene{j,7})
                    CondenTable(n,1:17) = {knownGene(j,1) GID(j) knownGene{j,2:7} ...
                        Location(k) Orientation(o) 'Intragenic' conden.Distance(a(m),1) knownGene{j,5}*conden.Distance(a(m),3) ...
                        conden.Window(a(m),1) conden.Window(a(m),2) conden.Reads(a(m),1) conden.Enrichment(a(m),1)};
                else
                    CondenTable(n,1:17) = {knownGene(j,1) GID(j) knownGene{j,2:7} ...
                        Location(k) Orientation(o) 'Intergenic' conden.Distance(a(m),1) knownGene{j,5}*conden.Distance(a(m),3) ...
                        conden.Window(a(m),1) conden.Window(a(m),2) conden.Reads(a(m),1) conden.Enrichment(a(m),1)};
                end
                fprintf(fid,sprintf('%s\t %d\t %s\t %s\t %s\t %d\t %d\t %d\t %s\t %s\t %s\t %d\t %d\t %d\t %d\t %d\t %d\n', ...
                    char(CondenTable{n,1}),cell2mat(CondenTable(n,2)),char(CondenTable{n,3}),char(CondenTable{n,4}), ...
                    char(CondenTable{n,5}),cell2mat(CondenTable(n,6)),cell2mat(CondenTable(n,7)), ...
                    cell2mat(CondenTable(n,8)),char(CondenTable{n,9}),char(CondenTable{n,10}), ...
                    char(CondenTable{n,11}),cell2mat(CondenTable(n,12)),cell2mat(CondenTable(n,13)), ...
                    cell2mat(CondenTable(n,14)),cell2mat(CondenTable(n,15)),cell2mat(CondenTable(n,16)),cell2mat(CondenTable(n,17))));
                n = n + 1;
            end
        end
    end
    if n < length(knownGene)
        CondenTable = CondenTable(1:n-1,:);
    end
    fclose(fid);clear fid
    
    assignin('base',sprintf('%s_%s_SigTable',headers{i},comment),SigTable);
    assignin('base',sprintf('%s_%s_CondenTable',headers{i},comment),CondenTable);
end

clear SigTable CondenTable n o a b m j i k Substruct l Orientation Columns ...
    GID Location GeneID C UCSC Win 