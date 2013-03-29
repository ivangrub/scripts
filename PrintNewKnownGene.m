clear,clc

fid = fopen('NewKnownGene.txt','w');

GeneID = fopen('GeneID','r');
known = fopen('../kgXref.txt','r');

k2r = textscan(known,'%s%s%s%s%s%s%s%s','delimiter','\t');
GI = textscan(GeneID,'%s%d','HeaderLines',1,'delimiter','\t');

for i = 1:length(k2r{1,1})
    a = find(strcmp(GI{1,1},k2r{1,1}(i)));
    if ~isempty(a)
        fprintf(fid,'%s\t%s\t%d\t%s\n',char(k2r{1,1}(i)),char(k2r{1,2}(i)),GI{1,2}(a),char(k2r{1,5}(i)));
    end
end
fclose(known);fclose(fid);fclose(GeneID);