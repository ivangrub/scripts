clear

person = {'Teppei' 'TY'};
headers = {'MEF_TBP' 'MEF_TAF7' 'MEF_PolII'};

for i = 1:length(headers)
    fid = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s_SPP_NarrowPeaks.txt',person{1},person{2},headers{i}),'r');
    fid2 = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s_SPP_NarrowPeaks.bed',person{1},person{2},headers{i}),'w');
    
    A = textscan(fid,'%s%d%d%s%d%s%f%f%f%f');
    fprintf(fid2,sprintf('track name=%s\n',headers{i}));
    for j = 1:length(A{1,1})
        fprintf(fid2,sprintf('%s\t%d\t%d\n',char(A{1,1}(j)),A{1,2}(j),A{1,3}(j)));
    end
    fclose(fid2);fclose(fid);clear fid fid2
end