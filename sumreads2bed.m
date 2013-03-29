%Print SumReads.mat file to a bed file

data = importdata('mm9.chr.length');
len = data.data;
chr = data.textdata;

for i = 1:length(headers)
    fprintf('Printing %s\n',headers{i});
    fid = fopen(sprintf('%s_%s_PostSPP.mm9.bed',perID{1},headers{i}),'w');
    x = eval(headers{i});
    for j = 1:length(chr)
        s = find(x.win.(chr{j})(2,:) ~= 0);
        for k = 1:length(s)
            if x.win.(chr{j})(1,s(k))+window-1 > len(j)
                yy = len(j);
            else
                yy = x.win.(chr{j})(1,s(k))+window-1;
            end
            fprintf(fid,'%s %d %d %s %d\n',chr{j},x.win.(chr{j})(1,s(k)),yy,headers{i},x.win.(chr{j})(2,s(k)));
        end
    end  
    fclose(fid);clear fid
end
    