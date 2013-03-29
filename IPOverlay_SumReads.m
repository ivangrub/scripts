%A scrheaderst to identify overlap between headers

clc

headers = {'Rad23b_Original' 'Oct4' 'Sox2'};
type = 'Lin';
mm9 = importdata('mm9.chr.length');

%% Identify which headerss you want to compare

chr = mm9.textdata;

for i = 1:length(headers)
    x = importdata(sprintf('%s_%s_Conden.txt',headers{i},type));
    assignin('base',sprintf('%s_IP_data',headers{i}),x.data);
    assignin('base',sprintf('%s_IP_text',headers{i}),x.textdata);
end

%% Set thresholds

Threshold = 500;                       %in units of bp

%% Calculate the overlap for the different combinations

m = 0;
for j = 1:length(headers)-1
    m =  m + length(headers) - j;
end

idx = 4:5;

TotalWin = zeros(1,m);
Overlap = zeros(1,m);

l = 1;
for i = 1:length(headers)-1
    for k = i+1:length(headers)
        total = 0;
        x = eval(sprintf('%s_IP_data',headers{i}));
        xx = eval(sprintf('%s_IP_text',headers{i}));
        y = eval(sprintf('%s_IP_data',headers{k}));
        yy = eval(sprintf('%s_IP_text',headers{k}));
        
        for j = 1:length(chr)
       
            XX = strcmp(chr{j},xx(2:end,5));
            YY = strcmp(chr{j},yy(2:end,5));
            
            A = unique(x(XX,idx),'rows');
            B = unique(y(YY,idx),'rows');
            total = total + length(A);
            
            it = 1;
            it2 = 1;
            for n = 1:length(A)
                Diff =  abs(A(n,1) - B(:,1)) <= Threshold;
                points = find(Diff == 1);
                if ~isempty(points)
                    OL2(j,it2:it2+length(points)-1) = points';
                    OL(j,it) = n;
                    it2 = it2 + length(points);
                    it = it + 1;
                end
            end
        end
        assignin('base',sprintf('%s_%s',headers{i},headers{k}),OL);
        assignin('base',sprintf('%s_%s_to2ndIP',headers{i},headers{k}),OL2);
        TotalWin(l) = total;
        Overlap(l) = length(find(OL));
        l = l + 1;
    end
end

Ratios = Overlap./TotalWin;

clear x y over o m l k j i