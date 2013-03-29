%Create a fasta file with the sequences from the called peaks for each IP. This will be used as an input
%for MUSCLE.
clear,clc

person = {'Teppei' 'TY'};
headers = {'ES_TBP' 'ES_PolII' 'ES_TAF1' 'ES_TAF7' 'MEF_TBP' 'MEF_PolII' 'MEF_TAF7'};

load knownGene.mm9.mat

data = importdata('mm9.chr.length');
chr = data.textdata;
len = data.data;

for k = 1:length(headers)
    fprintf('Printing %s\n',headers{k})
    fid = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s_SPP_PointPeaks.txt',person{1},person{2},headers{k}),'r');
    upstrea = fopen(sprintf('../ChIP-Seq Data/%s/FromPointPeaks_SPP_%s_Upstream_PeakSequence.fa',person{1},headers{k}),'w');
    tss = fopen(sprintf('../ChIP-Seq Data/%s/FromPointPeaks_SPP_%s_TSS_PeakSequence.fa',person{1},headers{k}),'w');
    body = fopen(sprintf('../ChIP-Seq Data/%s/FromPointPeaks_SPP_%s_GeneBody_PeakSequence.fa',person{1},headers{k}),'w');
    tes = fopen(sprintf('../ChIP-Seq Data/%s/FromPointPeaks_SPP_%s_TES_PeakSequence.fa',person{1},headers{k}),'w');
    
    A = textscan(fid,'%s%d%f%f%f%f%f','HeaderLines',1,'delimiter','\t');
    fclose(fid);clear fid
    for i = 1:length(chr)
        ind = find(strcmp(A{1,1}(:),chr{i}));
        c = find(strcmp(knownGene(:,4),chr{i}));
        Plus = find(strcmp(knownGene(:,4),chr{i}) & cell2mat(knownGene(:,5)) == 1);
        Minus = find(strcmp(knownGene(:,4),chr{i}) & cell2mat(knownGene(:,5)) == -1);
        x = fastaread(sprintf('MM9CHR/%s.fa',chr{i}));
        
        Peaks = cell2mat(A(1,2));
        B = double(Peaks(ind));
        fprintf('Printing %s\n',chr{i})
        for j = 1:length(B)
            
            % TSS
            p = (abs(B(j) - cell2mat(knownGene(Plus,6))) <= 500);  
            m = (abs(B(j) - cell2mat(knownGene(Minus,7))) <= 500);
            if sum(p) > 0 || sum(m) > 0
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
            p = (B(j) - cell2mat(knownGene(Plus,6)) >= -10000 & B(j) - cell2mat(knownGene(Plus,6)) < -500); 
            m = (B(j) - cell2mat(knownGene(Minus,7)) <= 10000 & B(j) - cell2mat(knownGene(Minus,7)) > 500);
            if sum(p) > 0 || sum(m) > 0
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
            p = B(j) >= cell2mat(knownGene(c,6)) & B(j) <= cell2mat(knownGene(c,7));
            if sum(p) > 0
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
            p = (B(j) - cell2mat(knownGene(Plus,6)) > 0 & B(j) - cell2mat(knownGene(Plus,6)) <= 10000);
            m = (B(j) - cell2mat(knownGene(Minus,7)) >= -10000 & B(j) - cell2mat(knownGene(Minus,7)) < 0);
            if sum(p) > 0 || sum(m) > 0
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
