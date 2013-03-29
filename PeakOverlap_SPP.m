%Peak Overlap from the PreImmune Text Files

clc

load knownGene.mm9.mat

perID = {'Teppei' 'TY'};
headers = {'ES_TAF7' 'ES_TBP' 'ES_PolII' 'ES_TAF1' 'MEF_TAF7' 'MEF_TBP' 'MEF_PolII' };
width = 500;                  %How far can peak windows be apart for it to be considered an overlap?
headersind = [1 2];


merging = 1;

Len = zeros(1,length(headersind));
for i = 1
    fid = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s_SPP_NarrowPeaks.txt',perID{1},perID{2},headers{headersind(i)}),'r');
    IP1 = textscan(fid,'%s%d%d%s%u%s%f%d%f%d');
    fclose(fid);clear fid
    
    Len(1,i) = length(IP1{1,1});
    chrsubset = cell(1,length(IP1{1,1}));
    subset = zeros(length(IP1{1,1}),4);
    
    for j = i+1:length(headersind)
        fid = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s_SPP_NarrowPeaks.txt',perID{1},perID{2},headers{headersind(i)}),'r');
        IP2 = textscan(fid,'%s%d%d%s%u%s%f%d%f%d');
        fclose(fid);clear fid
        
        Len(1,j) = length(IP2{1,1});
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
                
                
                subset(IND:IND+sum(Overlap)-1,[3 4]) = [IP2{1,2}(chr&win) IP2{1,3}(chr&win)];
                
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
            IP1{1,4} = IP1{1,3} - IP1{1,2};
            
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
        clear win subset chrsubset summitsubset cover RUN jj IP2
    end
    assignin('base',sprintf('%s_relativeto_%s',strcat(headers{headersind(2:end)}),headers{headersind(1)}),IP1);
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

clear IP1 IP2 LeftEdge RightEdge left right cover