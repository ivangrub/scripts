%Print SumReads.mat file to a bed file

data = importdata('dm3_Benjamin.chr.length');
len = data.data;
chr = data.textdata;

for i = 1:length(headers)
    fprintf('Printing %s\n',headers{i});
    x = eval(headers{i});
    for j = 1:length(chr)
<<<<<<< .mine
        fid = fopen(sprintf('%s_%s.%s.dm3_Benjamin.wig',perID{1},headers{i},chr{j}),'w');
=======
        fid = fopen(sprintf('%s_%s.%s.doublecheck.dm3.wig',perID{1},headers{i},chr{j}),'w');
>>>>>>> .r188
        fprintf(fid,sprintf('track type=wiggle_0 name="%s_%s" description="%s"\nvariableStep chrom=%s span=%d\n',headers{i},chr{j},headers{i},chr{j},window));
        s = find(x.win.(chr{j})(2,:) ~= 0);
        for k = 1:length(s)
            fprintf(fid,'%d %d\n',x.win.(chr{j})(1,s(k)),x.win.(chr{j})(2,s(k)));
        end
        fclose(fid);clear fid
    end 
    tar(sprintf('%s_%s.tgz',perID{1},headers{i}),'*.wig')
    delete(sprintf('%s_%s.*.dm3_Benjamin.wig',perID{1},headers{i})) 
end

