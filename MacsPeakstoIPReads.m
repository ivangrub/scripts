%Output IP read values for another IPs Macs peaks
clc
combination ='';
person = 'Claudia';

HeaderString = '%s\t%s\t%s';
String = '%s\t%d\t%d';
for i = 1:length(headers)
    HeaderString = strcat(HeaderString,'\t%s');
    String = strcat(String,'\t%d');
end
HeaderString = strcat(HeaderString,'\n');
String = strcat(String,'\n');

filename = 'Marson_IP_OSNT_PeakOverlap.txt';

fid = fopen(filename);
data = textscan(fid,'%s%d%d%s');
fclose(fid);
clear fid

OS = strcmp(data{1,4},sprintf('%s',combination)) | strcmp(data{1,4},sprintf('N%s',combination)) ...
    | strcmp(data{1,4},sprintf('N%sT',combination)) | strcmp(data{1,4},sprintf('%sT',combination));

Table = zeros(sum(OS),length(headers));
Chr = cell(1,sum(OS));
Win = zeros(sum(OS),2);
for j = 1:length(headers)
    x = eval(headers{j});
    k = 1;
    for i = 1:length(data{1,1})
        if ~OS(i)
            continue
        end
        chr = data{1,1}(i);
        window =  [data{1,2}(i) data{1,3}(i)];
        Table(k,j) = sum(x.bp.(chr{1}) >= window(1) & x.bp.(chr{1}) <= window(2));
        Chr(k) = chr;
        Win(k,1:2) = window;
        k = k + 1;
    end
end
if strcmp(combination,'')
    combination = 'NaN';
end
fid = fopen(sprintf('%s_Yick_w%s.txt',combination,person),'w');
fprintf(fid,HeaderString,'Chr','Start','End',headers{:});
for i = 1:length(Table)
    fprintf(fid,sprintf(String,Chr{i},Win(i,1),Win(i,2),Table(i,:)));
end
fclose(fid);

clear filename OS x j i chr window ind IPwindow data fid ans k Win Chr