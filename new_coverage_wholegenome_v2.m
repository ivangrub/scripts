%% Calculate the Coverage across the whole genome

%Run this script when you approximately know the length of the full genome
%analysis. New_coverage_wholegenome.m will be easier for first time use
%since it just conacates the arrays, but it does take longer.

clear,clc

person = {'Francisco','ESC'};
headers = {'TBP' 'TAF1' 'PolII_Ser5' 'PolII_N' 'PI' 'MCPI'};

fprintf('Opening the data files\n')
for i = 1:length(headers)
    x = open(sprintf('%s/%s_%s.mm9.mat',person{1},person{2},headers{i}));
    assignin('base',sprintf('data%d',i),x);
end

data = open('mm9length.mat');
mm9length = data.mm9length;

clear x i data 

%% Initialize the necessary matrices
window = 125;                               %In basepairs

fprintf('Initializing the matrices\n')
win = floor(window/25);
Chr = chr_save(win,data1,mm9length);

m = length(Chr);
n = length(headers);
Reads = zeros(m,n);
ReadsF = zeros(m,n);
ReadsR = zeros(m,n);
NewReads = zeros(1,n);
NewReadsF = zeros(1,n);
NewReadsR = zeros(1,n);
A = zeros(m,2);

%% Run the sliding window
yy = strrep(fieldnames(data1.chip),'chr','');
chr = cellstr(yy);

fprintf('Running the sliding window\n')
y = 1;
k = 1;
for i = 1:length(chr);
    fprintf('On chr%s\n',chr{i})
    x = floor(mm9length(i)/25);
    tss1 = 1:win:x-win;
    tss2 = win+1:win:x;
    a = length(tss1);
    C = [tss1' tss2'];
    A(y:y+a-1,1:2) = C;
    for j = 1:a
        for z = 1:n
            R = eval(sprintf('data%d',z));
            NewReads(z) = sum(R.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
            NewReadsF(z) = sum(R.chipF.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
            NewReadsR(z) = sum(R.chip.(sprintf('chr%s',chr{i}))(C(j,1):C(j,2)));
        end
        Reads(k,:) = NewReads;
        ReadsF(k,:) = NewReadsF;
        ReadsR(k,:) = NewReadsR;
        k = k + 1;
    end
    y = y + a;
end

FullGenome_Windows = A;
chr_num = yy;
%% Assemble the knownGenome annotations

fprintf('Saving knownGene\n')
data = open('knownGene_mm9.mat');
CpG = data.CpG;
knownGene = data.knownGene;

clear A C y tss1 tss2  data* yy NewReads a k chr i j x m NewReadsF NewReadsR ...
    person data n win window
%% Run a script to calculate the distances to the closest TSS

fprintf('Calculating the distance to the TSS\n')
Dist2TSS;

clear Chr