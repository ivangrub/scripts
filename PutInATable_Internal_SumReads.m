%Put the condensed/significant genes into a table

clc

%% Import the GeneID and set controls/comments

A = importdata('GeneID');
GeneID = A.data;
UCSC = A.textdata(2:end,1);
person = perID{1};

comment = 'Poisson';

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

%% Initialize the necessary matrices and options

Location = {'TSS' 'TSS' 'TES' 'TES'};
Orientation = {'Upstream' 'Downstream'};
Columns = [3 5 8 10];
Win = [1 1 6 6];
SigTable = cell(length(knownGene),17);
CondenTable = cell(length(knownGene),17);

%% Run through the significant and condensed structures


for i = 1:length(headers)
    l = 1;
    n = 1;
    
    sig = eval(sprintf('%s_Poiss_Sig',headers{i}));
    
    %Printing table for significant subset
    fprintf('Printing table for %s significant\n',headers{i})
    fid = fopen(sprintf('%s_%s_%s_Sig.txt',person,headers{i},comment),'w');
    fprintf(fid,'UCSC-ID\tGENE-ID\tGene Abbreviation\tGene Descript.\tChr\tStrand\tTSS\tTES\tFeature\tRelative Position\tInter or Intra\tMiddle of Window\tDistance to Feature\tLeft Edge\tRight Edge\tReads\tEnrichment\n');
    for j = 1:length(knownGene)
        [a,b] = find((j == sig.Distance(:,[2 4 7 9])));
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
                SigTable(l,1:16) = {knownGene(j,1) GID(j) knownGene{j,2:7} ...
                    Location(b(k)) Orientation(o) 'Intragenic' sig.Distance(a(k),Win(b(k))) knownGene{j,5}*sig.Distance(a(k),Columns(b(k))) ...
                    sig.Window(1,a(k)) sig.Window(2,a(k)) sig.Reads(a(k))};
            else
                SigTable(l,1:16) = {knownGene(j,1) GID(j) knownGene{j,2:7} ...
                    Location(b(k)) Orientation(o) 'Intergenic' sig.Distance(a(k),Win(b(k))) knownGene{j,5}*sig.Distance(a(k),Columns(b(k))) ...
                    sig.Window(1,a(k)) sig.Window(2,a(k)) sig.Reads(1,a(k))};
            end
            fprintf(fid,sprintf('%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n', ...
                char(SigTable{l,1}),cell2mat(SigTable(l,2)),char(SigTable{l,3}),char(SigTable{l,4}), ...
                char(SigTable{l,5}),cell2mat(SigTable(l,6)),cell2mat(SigTable(l,7)), ...
                cell2mat(SigTable(l,8)),char(SigTable{l,9}),char(SigTable{l,10}), ...
                char(SigTable{l,11}),cell2mat(SigTable(l,12)),cell2mat(SigTable(l,13)), ...
                cell2mat(SigTable(l,14)),cell2mat(SigTable(l,15)),cell2mat(SigTable(l,16))));
            l = l + 1;
        end
    end
    if l < length(knownGene)
        SigTable = SigTable(1:l-1,:);
    end
    fclose(fid);clear fid
    
    fprintf('Printing table for %s condensed\n',headers{i})
    %Printing table for condensed windows
    fid = fopen(sprintf('%s_%s_%s_Conden.txt',person,headers{i},comment),'w');
    fprintf(fid,'UCSC-ID\tGENE-ID\tGene Abbreviation\tGene Descript.\tChr\tStrand\tTSS\tTES\tRelative to TSS\tRelative to TES\tInter or Intra\tMiddle of Window\tTSS Dist\tTES Dist\tLeft Edge\tRight Edge\tReads\tEnrichment\n');
    conden = eval(sprintf('%s_Poiss_Conden',headers{i}));
    for j = 1:length(conden)
        if sign(conden(j,5)) == 1
            o = 2;
        else
            o = 1;
        end
        if sign(conden(j,6)) == 1
            l = 2;
        else
            l = 1;
        end
        if (conden(j,3) > knownGene{conden(j,4),6}) && ...
                (conden(j,3) < knownGene{conden(j,4),7})
            CondenTable(n,1:17) = {knownGene(conden(j,4),1) GID(conden(j,4)) knownGene{conden(j,4),2:7} ...
                Orientation(o) Orientation(l) 'Intragenic' conden(j,3) conden(j,5) conden(j,6)...
                conden(j,1) conden(j,2) conden(j,7)};
        else
            CondenTable(n,1:17) = {knownGene(conden(j,4),1) GID(conden(j,4)) knownGene{conden(j,4),2:7}...
                Orientation(o) Orientation(l) 'Intergenic' conden(j,3) conden(j,5) conden(j,6) ...
                conden(j,1) conden(j,2) conden(j,7)};
        end
        fprintf(fid,sprintf('%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n', ...
            char(CondenTable{n,1}),cell2mat(CondenTable(n,2)),char(CondenTable{n,3}),char(CondenTable{n,4}), ...
            char(CondenTable{n,5}),cell2mat(CondenTable(n,6)),cell2mat(CondenTable(n,7)), ...
            cell2mat(CondenTable(n,8)),char(CondenTable{n,9}),char(CondenTable{n,10}),char(CondenTable{n,11}), ...
            cell2mat(CondenTable(n,12)),cell2mat(CondenTable(n,13)),cell2mat(CondenTable(n,14)), ...
            cell2mat(CondenTable(n,15)),cell2mat(CondenTable(n,16)),cell2mat(CondenTable(n,17))));
        n = n + 1;
    end

if n < length(knownGene)
    CondenTable = CondenTable(1:n-1,:);
end
fclose(fid);clear fid

assignin('base',sprintf('%s_%s_Sig',headers{i},comment),SigTable);
assignin('base',sprintf('%s_%s_Conden',headers{i},comment),CondenTable);

end
clear SigTable CondenTable n o a b m j i k Substruct l Orientation Columns ...
    GID Location GeneID C UCSC Win 
