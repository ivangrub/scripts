%Peak Overlap from the PreImmune Text Files

clear,clc

load('knownGene.mm9.mat')
chr = importdata('mm9.chr.length');
Fields = chr.textdata;
clear chr

perID = {'TY'};
headers = {'TAF7' 'TBP' 'PolII_Ser5' 'PolII_N' 'PI' 'MEF_TAF7' 'MEF_TBP' 'MEF_PolII' 'MEF_PI'};
width = 500;        %How far can peak windows be apart for it to be considered an overlap?
headersind = [6 7 8];
controlind = 9;

merging = 1;

Len = zeros(1,length(headersind));
for i = 1
    A = importdata(sprintf('%s_%s_relativeto_%s_Lin_Conden.txt',perID{1},headers{headersind(i)},headers{controlind(1)}),'\t',1);
    data = A.data;
    text = A.textdata;
    IND = 1;
    for k = 1:length(Fields)
        chr = strcmp(Fields(k),text(2:end,5));
        correct = unique([data(chr,4) data(chr,5)],'rows');
        IP1{1,1}(IND:IND+length(correct)-1) = Fields(k);
        IP1{1,2}(IND:IND+length(correct)-1) = correct(:,1);
        IP1{1,3}(IND:IND+length(correct)-1) = correct(:,2);
        IND = IND + length(correct);
    end
    Len(1,i) = length(IP1{1,1});
    clear text data A
        
    chrsubset = cell(1,length(IP1{1,1}));
    subset = zeros(length(IP1{1,1}),4);
        
    for j = i + 1:length(headersind)
        A = importdata(sprintf('%s_%s_relativeto_%s_Lin_Conden.txt',perID{1},headers{headersind(j)},headers{controlind(1)}),'\t',1);
        data = A.data;
        text = A.textdata;
        IND = 1;
        for k = 1:length(Fields)
            chr = strcmp(Fields(k),text(2:end,5));
            correct = unique([data(chr,4) data(chr,5)],'rows');
            IP2{1,1}(IND:IND+length(correct)-1) = Fields(k);
            IP2{1,2}(IND:IND+length(correct)-1) = correct(:,1);
            IP2{1,3}(IND:IND+length(correct)-1) = correct(:,2);
            IND = IND + length(correct);
        end
        Len(1,j) = length(IP2{1,1});
        clear text data A
        
        IND = 1;
        for k = 1:length(IP1{1,1})
            chr = strcmp(IP2{1,1},IP1{1,1}(k));
            win = (IP2{1,2} >= IP1{1,2}(k)-width & IP2{1,3} <= IP1{1,3}(k) + width) | ...       %IP2 window fits inside of IP1 window
                (IP2{1,2} < IP1{1,2}(k)-width & IP2{1,3} >= IP1{1,2}(k) - width) | ...         %IP2 window straddles start of IP1 window
                (IP2{1,2} <= IP1{1,3}(k)+width & IP2{1,3} > IP1{1,3}(k) + width) | ...         %IP2 window staddles the end of IP1 window
                (IP2{1,2} <= IP1{1,2}(k)-width & IP2{1,3} >= IP1{1,3}(k) + width);              %IP1 window fits inside of the IP2 window
            Overlap = (chr & win);
            if sum(Overlap) == 1
                chrsubset(IND) = IP1{1,1}(k);
                subset(IND,[1 2 3 4]) = [IP1{1,2}(k) IP1{1,3}(k) IP2{1,2}(chr & win) IP2{1,3}(chr & win)];
                IND = IND + 1;
                continue
            elseif sum(Overlap) > 1
                
                chrsubset(IND:IND+sum(Overlap)-1) = IP1{1,1}(k);
                
                subset(IND:IND+sum(Overlap)-1,1) = IP1{1,2}(k);
                subset(IND:IND+sum(Overlap)-1,2) = IP1{1,3}(k);
                
                
                subset(IND:IND+sum(Overlap)-1,[3 4]) = [IP2{1,2}(chr&win)' IP2{1,3}(chr&win)'];
                IND = IND + sum(Overlap);
                
                continue
            elseif sum(Overlap) == 0
                chrsubset(IND) = IP1{1,1}(k);
                subset(IND,[1 2 3 4]) = [IP1{1,2}(k) IP1{1,3}(k) 0 0];
                IND = IND + 1;
                continue
            end
        end
        cover = subset(:,3) ~= 0;
        subset = subset(cover,:);
        win = zeros(sum(cover),4);
        clear IP1
        IP1{1,2} = zeros(sum(cover),1);
        IP1{1,3} = zeros(sum(cover),1);
        if merging == 1
            IP1{1,1} = chrsubset(cover);
            IP1{1,2} = min(subset,[],2);
            IP1{1,3} = max(subset,[],2);
        else
            IP1{1,1} = chrsubset(cover);
            win(:,1) = (subset(:,3) >= subset(:,1) & subset(:,4) <= subset(:,2)); 
            win(:,2) = (subset(:,3) < subset(:,1) & subset(:,4) >= subset(:,1));
            win(:,3) = (subset(:,3) <= subset(:,2) & subset(:,4) > subset(:,2)); 
            win(:,4) = (subset(:,3) <= subset(:,1) & subset(:,4) >= subset(:,2));
            win = logical(win);
            IP1{1,2}(win(:,1)) = subset(win(:,1),3);
            IP1{1,3}(win(:,1)) = subset(win(:,1),4);
            IP1{1,2}(win(:,2)) = subset(win(:,2),1);
            IP1{1,3}(win(:,2)) = subset(win(:,2),4);
            IP1{1,2}(win(:,3)) = subset(win(:,3),3);
            IP1{1,3}(win(:,3)) = subset(win(:,3),2);
            IP1{1,2}(win(:,4)) = subset(win(:,4),1); 
            IP1{1,3}(win(:,4)) = subset(win(:,4),2);
        end
        clear win subset chrsubset summitsubset cover IP2
    end
    assignin('base',sprintf('%s_relativeto_%s',strcat(headers{headersind(2:end)}),headers{headersind(i)}),IP1);
end

subset = [IP1{1,2} IP1{1,3}];
chrsubset = IP1{1,1};

[m,n] = size(subset);

InGene = zeros(1,m);
for i = 1:m
    knownGenesubset = strcmp(chrsubset(i),knownGene(:,4));
    left = subset(i,1);
    right = subset(i,2);
    LeftEdge = cell2mat(knownGene(knownGenesubset,6))-1000;
    RightEdge = cell2mat(knownGene(knownGenesubset,7))+1000;
    win = (left >= LeftEdge & right <= RightEdge) | ...       %IP2 window fits inside of OS window
        (left < LeftEdge & right >= LeftEdge) | ...         %IP2 window straddles start of OS window
        (left <= RightEdge & right > RightEdge) | ...         %IP2 window staddles the end of the window
        (left <= LeftEdge & right >= RightEdge);              %OS window fits inside of the IP2 window
    InGene(1,i) = sum(win);
end
Total = sum(InGene);