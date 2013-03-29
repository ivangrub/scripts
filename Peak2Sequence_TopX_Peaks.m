%Create a fasta file with the sequences from the called peaks for each IP. This will be used as an input 
%for MUSCLE.
clear, clc

Top = 100;

person = {'Teppei' 'TY'};
headers = {'ES_TBP' 'ES_PolII' 'ES_TAF1' 'ES_TAF7' 'MEF_TBP' 'MEF_PolII' 'MEF_TAF7'};

data = importdata('mm9.chr.length');
len = data.data;
chr = data.textdata;

for k = 1:length(headers)
    fid = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s_SPP_NarrowPeaks.txt',person{1},person{2},headers{k}),'r');  
    fid2 = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s_Top%d_SPP_PeakSequence.fa',person{1},person{2},headers{k},Top),'w');
    A = textscan(fid,'%s%d%d%s%d%s%f%f%f%f');
    
    for i = 1:length(chr)
        B = find(strcmp(A{1,1}(1:Top),chr{i}));
        x = fastaread(sprintf('~/Documents/MATLAB/Tjian Lab/ChIP-Seq/MM9CHR/%s.fa',chr{i}));
        mid = round(mean([A{1,2}(B) A{1,3}(B)],2));
        
        fprintf('Printing %s\n',chr{i})
        for j = 1:length(mid)
            if mid(j) > 125 && mid(j)+125 < len(i)
                position = sprintf('%s!%d!%d',char(chr{i}),mid(j)-125,mid(j)+125);
                fprintf(fid2,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(mid(j)-125:mid(j)+125))));
            elseif mid(j) < 125
                position = sprintf('%s!%d!%d',char(chr{i}),mid(j),mid(j)+250);
                fprintf(fid2,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(mid(j):mid(j)+250))));
            elseif mid(j)+125 > len(i)
                position = sprintf('%s!%d!%d',char(chr{i}),mid(j)-250,len(i));
                fprintf(fid2,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(mid(j)-250:len(i)))));
            end
        end
    end
    fclose(fid);fclose(fid2);clear fid fid2
end
