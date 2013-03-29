clear,clc

fid = fopen('NewKnownGene.txt');
table = textscan(fid,'%s%s%d%s','delimiter','\t');
fclose(fid);clear fid;

fid = fopen('KnownGene_wTSS.txt','w');
[x,y,z] = unique(table{3});

load knownGene.mm9.mat

for i = 1:length(y)
    id = find(strcmp(table{4}(y(i)),knownGene(:,2)));
    if isempty(id)
       continue
    end
    if cell2mat(knownGene(id(1),5)) == -1
        strand = '-';
    else
        strand = '+';
    end
    for j = 1:length(id)
        fprintf(fid,'%s\t%s\t%d\t%s\t%s\t%s\t%d\t%d\n',char(table{1}(y(i))),char(table{2}(y(i))),table{3}(y(i)),char(table{4}(y(i))),char(knownGene{id(1),4}),strand,cell2mat(knownGene(id(j),6)),cell2mat(knownGene(id(j),7)));
    end
    
end
fclose(fid);clear fid