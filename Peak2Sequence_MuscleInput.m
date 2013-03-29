%Create a fasta file with the sequences from the called peaks for each IP. This will be used as an input 
%for MUSCLE.
clear, clc

person = {'James' 'James'};
headers = {'TAF3'};

data = importdata('mm9.chr.length');
len = data.data;
chr = data.textdata;

for k = 1:length(headers)
    fid = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s_SPP_PointPeaks.txt',person{1},person{2},headers{k}),'r');  
    fid2 = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s_SPP_PeakSequence.fa',person{1},person{2},headers{k}),'w');
    A = textscan(fid,'%s%d%f%f%f%f%f','HeaderLines',1,'delimiter','\t');
    
    for i = 1:length(chr)
        B = find(strcmp(A{1,1}(:),chr{i}));
        x = fastaread(sprintf('MM9CHR/%s.fa',chr{i}));
        
        
        fprintf('Printing %s\n',chr{i})
        for j = 1:length(B)
            if A{1,2}(B(j)) > 125 && A{1,2}(B(j))+125 < len(i)
                position = sprintf('%s!%d!%d',char(A{1,1}(B(j))),A{1,2}(B(j))-125,A{1,2}(B(j))+125);
                fprintf(fid2,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(A{1,2}(B(j))-125:A{1,2}(B(j))+125))));
            elseif A{1,2}(B(j)) < 125
                position = sprintf('%s!%d!%d',char(A{1,1}(B(j))),A{1,2}(B(j)),A{1,2}(B(j))+250);
                fprintf(fid2,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(A{1,2}(B(j)):A{1,2}(B(j))+250))));
            elseif A{1,2}(B(j))+125 > len(i)
                position = sprintf('%s!%d!%d',char(A{1,1}(B(j))),A{1,2}(B(j))-250,len(i));
                fprintf(fid2,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(A{1,2}(B(j))-250:len(i)))));
            end
        end
    end
end
