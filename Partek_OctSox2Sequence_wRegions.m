clear,clc

fid = fopen('Rad23b.1_NoOctSox_45064.txt');
fid2 = fopen('Rad23b.1_OctSox_14705.txt');

headers = {'NoOS' 'OS'};

OS = textscan(fid2,'%s%d%d','HeaderLines',1,'delimiter','\t');
NoOS = textscan(fid,'%s%d%d','HeaderLines',1,'delimiter','\t');

load knownGene.mm9.mat

fclose(fid);fclose(fid2);clear fid fid2

for i = 1:length(headers)
    fprintf('Printing %s\n',headers{i})
    x = eval(headers{i});
    mid = round(mean([x{1,2} x{1,3}],2));
    
    upstrea = fopen(sprintf('Partek_%s_Upstream_PeakSequence.fa',headers{i}),'w');
    tss = fopen(sprintf('Partek_%s_TSS_PeakSequence.fa',headers{i}),'w');
    body = fopen(sprintf('Partek_%s_GeneBody_PeakSequence.fa',headers{i}),'w');
    tes = fopen(sprintf('Partek_%s_TES_PeakSequence.fa',headers{i}),'w');
    
    for j = 1:length(mid)
        chr = strcat('chr',x{1,1}(j));
        
        if j == 1 || ~strcmp(x{1,1}(j),x{1,1}(j-1))
            y = fastaread(sprintf('MM9CHR/%s.fa',chr{1}));
            len = length(y.Sequence);
            fprintf('Printing %s\n',chr{1})
        end
        
        c = find(strcmp(knownGene(:,4),chr{1}));
        Plus = find(strcmp(knownGene(:,4),chr{1}) & cell2mat(knownGene(:,5)) == 1);
        Minus = find(strcmp(knownGene(:,4),chr{1}) & cell2mat(knownGene(:,5)) == -1);
        
        
        
        % TSS
        p = (abs(mid(j) - cell2mat(knownGene(Plus,6))) <= 500);
        m = (abs(mid(j) - cell2mat(knownGene(Minus,7))) <= 500);
        if sum(p) > 0 || sum(m) > 0
            if mid(j) > 125 && mid(j)+125 < len
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j)-125,mid(j)+125);
                fprintf(tss,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j)-125:mid(j)+125))));
            elseif mid(j) < 125
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j),mid(j)+250);
                fprintf(tss,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j):mid(j)+250))));
            elseif mid(j)+125 > len
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j)-250,len);
                fprintf(tss,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j)-250:len))));
            end
        end
        
        %Upstream
        p = (mid(j) - cell2mat(knownGene(Plus,6)) >= -10000 & mid(j) - cell2mat(knownGene(Plus,6)) < -500);
        m = (mid(j) - cell2mat(knownGene(Minus,7)) <= 10000 & mid(j) - cell2mat(knownGene(Minus,7)) > 500);
        if sum(p) > 0 || sum(m) > 0
            if mid(j) > 125 && mid(j)+125 < len
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j)-125,mid(j)+125);
                fprintf(upstrea,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j)-125:mid(j)+125))));
            elseif mid(j) < 125
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j),mid(j)+250);
                fprintf(upstrea,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j):mid(j)+250))));
            elseif mid(j)+125 > len
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j)-250,len);
                fprintf(upstrea,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j)-250:len))));
            end
        end
        
        %Body
        p = mid(j) >= cell2mat(knownGene(c,6)) & mid(j) <= cell2mat(knownGene(c,7));
        if sum(p) > 0
            if mid(j) > 125 && mid(j)+125 < len
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j)-125,mid(j)+125);
                fprintf(body,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j)-125:mid(j)+125))));
            elseif mid(j) < 125
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j),mid(j)+250);
                fprintf(body,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j):mid(j)+250))));
            elseif mid(j)+125 > len
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j)-250,len);
                fprintf(body,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j)-250:len))));
            end
        end
        
        %TES
        p = (mid(j) - cell2mat(knownGene(Plus,6)) > 0 & mid(j) - cell2mat(knownGene(Plus,6)) <= 10000);
        m = (mid(j) - cell2mat(knownGene(Minus,7)) >= -10000 & mid(j) - cell2mat(knownGene(Minus,7)) < 0);
        if sum(p) > 0 || sum(m) > 0
            if mid(j) > 125 && mid(j)+125 < len
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j)-125,mid(j)+125);
                fprintf(tes,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j)-125:mid(j)+125))));
            elseif mid(j) < 125
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j),mid(j)+250);
                fprintf(tes,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j):mid(j)+250))));
            elseif mid(j)+125 > len
                position = sprintf('%s!%d!%d',char(chr{1}),mid(j)-250,len);
                fprintf(tes,sprintf('>%s\n%s\n\n',position,upper(y.Sequence(mid(j)-250:len))));
            end
        end
    end
    fclose(body);fclose(tss);fclose(tes);fclose(upstrea);clear body tss tes upstre
end


