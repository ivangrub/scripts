%Create a fasta file with the sequences from the called peaks for each IP. This will be used as an input
%for MUSCLE.
clear,clc

person = {'James' 'James'};
headers = {'TAF3'};

fid = fopen('known_wTSS.txt');
known = textscan(fid,'%s%s%d%s%s%s%d%d','delimiter','\t');
fclose(fid);clear fid

data = importdata('mm9.chr.length');
chr = data.textdata;
len = data.data;
sections = {'TSS' 'Upstream' 'Body' 'TES'};

for k = 1:length(headers)
    fprintf('Printing %s\n',headers{k})
    fid = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s_SPP_PointPeaks.txt',person{1},person{2},headers{k}),'r');
    upstrea = fopen(sprintf('../ChIP-Seq Data/%s/FromPointPeaks_SPP_%s_Upstream_PeakSequence.SingleGene.fa',person{1},headers{k}),'w');
    tss = fopen(sprintf('../ChIP-Seq Data/%s/FromPointPeaks_SPP_%s_TSS_PeakSequence.SingleGene.fa',person{1},headers{k}),'w');
    body = fopen(sprintf('../ChIP-Seq Data/%s/FromPointPeaks_SPP_%s_GeneBody_PeakSequence.SingleGene.fa',person{1},headers{k}),'w');
    tes = fopen(sprintf('../ChIP-Seq Data/%s/FromPointPeaks_SPP_%s_TES_PeakSequence.SingleGene.fa',person{1},headers{k}),'w');
    
    A = textscan(fid,'%s%d%f%f%f%f%f','HeaderLines',1,'delimiter','\t');
    fclose(fid);clear fid
    for i = 1:length(chr)
        ind = find(strcmp(A{1,1}(:),chr{i}));
        c = find(strcmp(known{5},chr{i}));
        Plus = find(strcmp(known{5},chr{i}) & strcmp(known{6},'+'));
        Minus = find(strcmp(known{5},chr{i}) & strcmp(known{6},'-'));
        x = fastaread(sprintf('MM9CHR/%s.fa',chr{i}));
        
        Peaks = cell2mat(A(1,2));
        B = double(Peaks(ind));
        fprintf('Printing %s\n',chr{i})
        for j = 1:length(B)
            conditions = zeros(4,1);
            promoter_plus = (abs(B(j) - cell2mat(known{6}(Plus))) <= 500);  
            promoter_minus = (abs(B(j) - cell2mat(known{8}(Minus))) <= 500);
            if sum([length(unique(known{4}(promoter_plus))) length(unique(known{4}(promoter_minus)))]) == 1
                conditions(1) = 1;
            end
            conditions = [sum((abs(B(j) - cell2mat(known{7}(Plus))) <= 500)) sum((abs(B(j) - cell2mat(known{7}(Minus))) <= 500));...
                sum((B(j) - cell2mat(known{7}(Plus)) >= -10000 & B(j) - cell2mat(known{7}(Plus)) < -500)) sum((B(j) - cell2mat(known{8}(Minus)) <= 10000 & B(j) - cell2mat(known{8}(Minus)) > 500));...
                sum(B(j) >= cell2mat(known{7}(c)) & B(j) <= cell2mat(known{8}(c))) 0;...
                sum((B(j) - cell2mat(known{7}(Plus)) > 0 & B(j) - cell2mat(known{7}(Plus)) <= 10000)) sum((B(j) - cell2mat(known{8}(Minus)) >= -10000 & B(j) - cell2mat(known{8}(Minus)) < 0))];
            
            if sum(sum(conditions)) == 1
                index = find(sum(conditions,2));
                sect = sections(index);
            else
                continue
            end
            % TSS
            
            if strcmp(sect,'TSS')
                if B(j) > 125 && B(j)+125 < len(i)
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j)-125,B(j)+125);
                    fprintf(tss,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j)-125:B(j)+125))));
                elseif B(j) < 125
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j),B(j)+250);
                    fprintf(tss,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j):B(j)+250))));
                elseif B(j)+125 > len(i)
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j)-250,len(i));
                    fprintf(tss,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j)-250:len(i)))));
                end
            end
            
            %Upstream
            
            if strcmp(sect,'Upstream')
                if B(j) > 125 && B(j)+125 < len(i)
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j)-125,B(j)+125);
                    fprintf(upstrea,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j)-125:B(j)+125))));
                elseif B(j) < 125
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j),B(j)+250);
                    fprintf(upstrea,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j):B(j)+250))));
                elseif B(j)+125 > len(i)
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j)-250,len(i));
                    fprintf(upstrea,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j)-250:len(i)))));
                end
            end
            
            %Body
            
            if strcmp(sect,'Body')
                if B(j) > 125 && B(j)+125 < len(i)
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j)-125,B(j)+125);
                    fprintf(body,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j)-125:B(j)+125))));
                elseif B(j) < 125
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j),B(j)+250);
                    fprintf(body,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j):B(j)+250))));
                elseif B(j)+125 > len(i)
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j)-250,len(i));
                    fprintf(body,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j)-250:len(i)))));
                end
            end
            
            %TES
            
            if strcmp(sect,'TES')
                if B(j) > 125 && B(j)+125 < len(i)
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j)-125,B(j)+125);
                    fprintf(tes,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j)-125:B(j)+125))));
                elseif B(j) < 125
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j),B(j)+250);
                    fprintf(tes,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j):B(j)+250))));
                elseif B(j)+125 > len(i)
                    position = sprintf('%s!%d!%d',char(chr{i}),B(j)-250,len(i));
                    fprintf(tes,sprintf('>%s\n%s\n\n',position,upper(x.Sequence(B(j)-250:len(i)))));
                end
            end
        end
    end
    fclose(body);fclose(tss);fclose(tes);fclose(upstrea);clear body tss tes upstre
end
