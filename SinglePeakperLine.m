
headers = {'TAF7L_Pre' 'TAF7L_Post' 'Pparg_Pre' 'Pparg_Post'};

for n = 1:length(headers)
    
    peaks = fopen(sprintf('~/Desktop/HY_%s.Tommy.txt',headers{n}));
    IP = textscan(peaks,'%s%s%d%d%f%f%s%s%f%s%s%s%d','HeaderLines',1,'delimiter','\t');
    fclose(peaks);clear peaks
    
    line = fopen(sprintf('~/Desktop/HY_%s_SinglePeakLines.txt',headers{n}),'w');
    for i = 1:length(IP{1,7})
        str = char(IP{1,7}(i));
        l = strfind(str,'X');
        
        k = 1;
        for j = 1:length(l);
            fprintf(line,'%s\t%d\t%d\t%s\n',char(IP{1,2}(i)),IP{1,3}(i),IP{1,4}(i),str(k:l(j)-1));
            k = l(j) + 1;
        end
    end
    fclose(line);
end