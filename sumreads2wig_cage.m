%Print SumReads.mat file to a wig file

perID = {'Lea'};
headers = {'Cage_MinData_MacroSpecific' 'H4Ac' 'H3'};
data = importdata('mm9.chr.length');
len = data.data;
chr = data.textdata;

for i = 1:length(headers)
    fprintf('Printing %s\n',headers{i});
    x = eval(headers{i});
    for j = 1:length(chr)
        vec = 1:10:len(j)-1;
        fid = fopen(sprintf('%s_%s_ForCage.%s.mm9.wig',perID{1},headers{i},chr{j}),'w');
        fprintf(fid,sprintf('track type=wiggle_0 name="%s_%s" description="%s"\nvariableStep chrom=%s span=10\n',headers{i},chr{j},headers{i},chr{j}));
        for k = 1:length(vec)-1
            s = sum(x.bp.(chr{j}) >= vec(k) & x.bp.(chr{j}) < vec(k+1));
            if s > 0
                fprintf(fid,'%d %d\n',vec(k),s);
            end
        end
        fclose(fid);clear fid   
    end  
end
    