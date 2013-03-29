%Calculate the Overlap of Narrow SPP peaks with relation to the first IP
%that is inputted. Everything is taken with regards to IP{1}(1). Enter 1 if
%you want them to overlap, 0 if it is independent or NaN if you dont care

clear,clc

perID = {{'express' 'PostExpress'} {'express' 'TY'}};
IP = {{'TBP'} {'MEF_TBP'}};
OVERLAPCODE = 1;
percentile = [1 .25];

k = 1;
for i = 1:length(perID)
    for j = 1:length(IP{i})
        fid = fopen(sprintf('~/ChIP-Data/%s/%s_%s_SPP_NarrowPeaks.txt',char(perID{i}(1)),char(perID{i}(2)),char(IP{i}(j))),'r');
        data = textscan(fid,'%s%d%d%s%u%s%f%d%f%d');
        data{1} = int8(char(data{1}));
        fclose(fid);clear fid
        assignin('base',char(IP{i}(j)),data);
        headers(k) = IP{i}(j);
        k = k + 1;
    end
end

firstIP = eval(headers{1});
for j = 1:length(firstIP)
    if percentile(1) == 1
        firstIP{j} = firstIP{j}(1:round(percentile(2)*length(firstIP{j})),:);
    else
        firstIP{j} = firstIP{j}(round(length(firstIP{j})*percentile(1)):round(length(firstIP{j})*percentile(2)),:);
    end
end
chrfirstIP = firstIP{1};
PossibleOverlaps = length(firstIP{1});
TotalOverlaps = zeros(1,PossibleOverlaps);
Dist = zeros(1,PossibleOverlaps);
Peaks = zeros(1,length(headers)-1);
Sum = Peaks;
total_freq = zeros(2*(length(headers)-2),200);

%compared = 1;
zz = 1;
for i = 2:length(headers)
    if strcmp(headers{i},headers{1})
        Peaks(i-1) = length(firstIP{1});
        continue
    end
%     if isnan(OVERLAPCODE(i-1))
%         continue
%     end
    saveind(i-1) = i;
    nextIP = eval(headers{i});
    for j = 1:length(nextIP)
        if percentile(1) == 1
            nextIP{j} = nextIP{j}(1:round(length(nextIP{j})*percentile(2)),:);
        else
            nextIP{j} = nextIP{j}(round(length(nextIP{j})*percentile(1)):round(length(nextIP{j})*percentile(2)),:);
        end
    end
    chrnextIP = nextIP{1};
    chr = zeros(length(chrnextIP),2);

    for k = 1:length(TotalOverlaps)
%         if compared >= 2 && TotalOverlaps(k) == 0
%             continue 
%         end
        for z = 4:5
            chr(:,z-3) = chrfirstIP(k,z) == chrnextIP(:,z);
        end
        chrind = sum(chr,2) == 2;
        win = (firstIP{2}(k) >= nextIP{2} & firstIP{3}(k) <= nextIP{3}) | ...       %IP window fits inside of OS window
            (firstIP{2}(k) < nextIP{2} & firstIP{3}(k) >= nextIP{2}) | ...          %IP window straddlMEFstart of OS window
            (firstIP{2}(k) <= nextIP{3} & firstIP{3}(k) > nextIP{3}) | ...          %IP window staddlMEFthe end of the window
            (firstIP{2}(k) <= nextIP{2} & firstIP{3}(k) >= nextIP{3});              %OS window fits inside of the IP window
        
        Overlap = (chrind & win);
        if OVERLAPCODE(1) == 0
            set = 0;
            opp = 1;
        else
            set = 1;
            opp = 0;
        end
        if sum(Overlap) >= 1
            TotalOverlaps(k) = set;
            s = find(Overlap);
            mid = mean([firstIP{2}(k) firstIP{3}(k)],2);
            d = zeros(1,length(s));
            for j = 1:length(s)
                d(j) = mid - mean([nextIP{2}(s(j)) nextIP{3}(s(j))],2);
            end
            Dist(k) = mean(d);
        else
            TotalOverlaps(k) = opp;
        end       
    end
    p = logical(TotalOverlaps);
    
    [freq,xout] = hist(Dist(p),200);
    
    f = freq/sum(TotalOverlaps);
    total_freq(zz,:) = xout;
    total_freq(zz+1,:) = f;
    
    Peaks(i-1) = sum(TotalOverlaps);
   	Sum(i-1) = sum(abs(Dist(p)) <= 100)/sum(TotalOverlaps);
   	zz = zz + 2;
    %compared = compared + 1;
end
[m,n] = size(total_freq);
color = {'k','r','b','g'};
j = 1;
for k = 1:2:m
    plot(total_freq(k,:),total_freq(k+1,:),color{j})
    hold on
    j = j + 1;
end
axis([-1000 1000 0 1.1*max(max(total_freq(2:2:end,:)))])
xlabel('Distance (bp)')
ylabel(sprintf('Fraction of Top %d%% Overlapping Peaks',percentile(2)*100))
hold off
saveas(gcf,sprintf('~/ChIP-Data/%s/SummitPlot_oftop%d_relativeto_%s.pdf',char(perID{1}(1)),percentile(2)*100,headers{1}),'pdf');
Peaks
PercentagePeaks = Peaks/PossibleOverlaps*100
Sum

fprintf('Print Overlapped Annotation\n')
A = [firstIP{2}(p) firstIP{3}(p)];
annotate_overlap(char(chrfirstIP(p,:)),A,char(perID{1}(1)),headers{1},percentile(2)*100,'yes')
fprintf('Print Non-overlapped Annotation\n')
A = [firstIP{2}(p == 0) firstIP{3}(p == 0)];
annotate_overlap(char(chrfirstIP(p == 0,:)),A,char(perID{1}(1)),headers{1},percentile(2)*100,'no')