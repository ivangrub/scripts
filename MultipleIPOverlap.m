%Compare overlaps of multiple IPs

clear,clc

IP = {'MEF_TBP' 'MEF_PolII' };
relto = {'MEF_TAF7'};
Overlap = [1 1];

for i = 1:length(IP)
    name = sprintf('OverlapPeaksFor_%s_relativeto_%s.txt',IP{i},relto{1});
    fid = fopen(name);
    data = textscan(fid,'%s%d%d%s%d%d%d%d%d%d','HeaderLines',1);
    assignin('base',IP{i},data);
    fclose(fid);clear fid
end

x = eval(IP{1});
if Overlap(1) == 1
    subset = ~strcmp(x{1,4},'NaN');
else
    subset = strcmp(x{1,4},'NaN');
end

chr = x{1,1}(subset);
start = double(x{1,2}(subset));
endit = double(x{1,3}(subset));
for i = 2:length(IP)
    y = eval(IP{i});
    ysub = zeros(length(y{1,1}),sum(subset));
    ind = 1;
    for k = 1:sum(subset)
        sub = strcmp(y{1,1},chr{k}) & (start(k) == y{1,2}) & (endit(k) == y{1,3});
        if sum(sub) > 1
            for j = 1:sum(sub)
                ysub(:,ind) = sub(j);
                ind = ind + 1;
            end
        else
            ysub(:,ind) = sub;
            ind = ind + 1;
        end
    end
    clear subset
    [~,n] = size(ysub);
    ind = 1;
    for k = 1:n
        if Overlap(i) == 1
            sub = ~strcmp(y{1,4}(logical(ysub(:,k))),'NaN');
        else
            sub = strcmp(y{1,4}(logical(ysub(:,k))),'NaN');
        end
        if sum(sub) > 1
            for j = 1:sum(sub)
                subset(:,ind) = sub(j);
                ind = ind + 1;
            end
        elseif sum(sub) == 1
            subset(:,ind) = sub;
            ind = ind + 1;
        else
            subset(:,ind) = 0;
            ind = ind + 1;
        end
    end
end


TotalOverlappingPeaks = sum(subset);
FractionPeaks = TotalOverlappingPeaks/length(y{1,4});
clear data name fid x