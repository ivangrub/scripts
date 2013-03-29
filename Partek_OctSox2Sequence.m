clear,clc

fid = fopen('Rad23b.1_NoOctSox_45064.txt');
fid2 = fopen('Rad23b.1_OctSox_14705.txt');

headers = {'NoOS' 'OS'};

OS = textscan(fid2,'%s%d%d','HeaderLines',1,'delimiter','\t');
NoOS = textscan(fid,'%s%d%d','HeaderLines',1,'delimiter','\t');

fclose(fid);fclose(fid2);clear fid fid2

fid1 = fopen('Partek_Rad23b_NoOctSox.PeakSequence.fa','w');
fid2 = fopen('Partek_Rad23b_WithOctSox.PeakSequence.fa','w');

for i = 1:length(headers)
    fprintf('Printing %s\n',headers{i})
    x = eval(headers{i});
    mid = round(mean([x{1,2} x{1,3}],2));
    
    for j = 1:length(mid)
        chr = strcat('chr',x{1,1}(j));
        
        if j == 1 || ~strcmp(x{1,1}(j),x{1,1}(j-1))
            y = fastaread(sprintf('MM9CHR/%s.fa',chr{1}));
            len = length(y.Sequence);
        end
        
        
        if mid(j) > 125 && mid(j)+125 < len
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j)-125,mid(j)+125);
                fprintf(eval(sprintf('fid%d',i)),sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j)-125:mid(j)+125))));
            elseif mid(j) < 125
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j),mid(j)+250);
                fprintf(eval(sprintf('fid%d',i)),sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j):mid(j)+250))));
            elseif mid(j)+125 > len
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j)-250,len);
                fprintf(eval(sprintf('fid%d',i)),sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j)-250:len))));
        end
    end
end
fclose(fid1);fclose(fid2);clear fid fid2