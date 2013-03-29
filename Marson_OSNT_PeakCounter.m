combo = {'OS' 'O' 'S' ''};

filename = 'Marson_IP_OSNT_PeakOverlap.txt';

fid = fopen(filename);
data = textscan(fid,'%s%d%d%s');
fclose(fid);
clear fid

Total = zeros(1,4);
for i = 1:length(combo)
    OS = strcmp(data{1,4},sprintf('%s',combo{i})) | strcmp(data{1,4},sprintf('N%s',combo{i})) ...
        | strcmp(data{1,4},sprintf('N%sT',combo{i})) | strcmp(data{1,4},sprintf('%sT',combo{i}));
    Total(i) = sum(OS);
end

clear combo filename data OS i ans