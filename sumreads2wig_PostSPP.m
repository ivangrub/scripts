%Print SumReads.mat file to a wig file

data = importdata('mm9.chr.length');
len = data.data;
chr = data.textdata;

for i = 1:length(headers)
    fprintf('Printing %s\n',headers{i});
    x = eval(headers{i});
    for j = 1:length(chr)
        fid = fopen(sprintf('%s_%s_PostSPP.%s.mm9.wig',perID{1},headers{i},chr{j}),'w');
        fprintf(fid,sprintf('track type=wiggle_0 name="%s_%s" description="%s"\nvariableStep chrom=%s span=250\n',headers{i},chr{j},headers{i},chr{j}));
        s = find(x.win.(chr{j})(2,:) ~= 0);
        for k = 1:length(s)
            fprintf(fid,'%d %d\n',x.win.(chr{j})(1,s(k)),x.win.(chr{j})(2,s(k)));
        end
        fclose(fid);clear fid   
    end  
end
    