%A script to identify overlap between IP

clc
%% Identify which IPs you want to compare

if iscell(IP) == 0
    error('Create an IP cell with the IPs of interest')
end

for i = 1:length(IP)
    x = eval(sprintf('%s_%s_CondenTable',IP{i},comment));
    assignin('base',IP{i},x);
end

%% Set thresholds

Threshold = 1000;                       %in units of bp

%% Calculate the overlap for the different combinations

m = 0;
for j = 1:length(IP)-1
    m =  m + length(IP) - j;
end

TotalWin = zeros(m,1);
Overlap = zeros(m,1);

l = 1;
for i = 1:length(IP)-1
    for k = i+1:length(IP)
        o = 0;
        x = eval(IP{i});
        A = unique(cell2mat(x(:,[12 15:16])),'rows');
        y = eval(IP{k});
        [YY,idx] = unique(cell2mat(y(:,[12 15:16])),'rows');
        
        OL = zeros(length(A),1);
        j = 1;
        it = 1;
        while j <= length(A)
            same_win1 = find(cell2mat(x(j,12)) == cell2mat(x(j:j+10,12)),1,'last');
            same_win1 = same_win1 + j;
            Diff =  abs(cell2mat(x(same_win1,12)) - YY(:,1)) <= Threshold;
            points = find(Diff == 1);
            if isempty(points) == 0
                if length(points) > 1
                   chrom = 0;
                   n = 1;
                   while chrom == 0
                       if n <= length(points)
                            chrom = strcmp(x{j,5},y{idx(points(n)),5});
                            n = n + 1;
                       else
                           chrom = 2;
                       end
                   end
                   if chrom == 1
                       OL(it,1) = points(n-1);
                       o = o + 1;
                   end
                else
                    if strcmp(x{j,5},y{idx(points),5})
                        OL(it,1) = points;
                        o = o + 1;
                    end
                end
            end
            it = it + 1;
            j = same_win1 + 1;
        end
        assignin('base',sprintf('%s_%s',IP{i},IP{k}),OL);
        TotalWin(l) = length(A);
        Overlap(l) = o;
        l = l + 1;
    end
end

Ratios = Overlap./TotalWin;
    
clear x y over o m l k j i Threshold