%% Full Genome Data Normalization and Analysis

clc
controlidx = 3; 
IPidx = 1;

%% Initialize matrices

conden_len = zeros(1,35);

%% Process Structured Input
x = eval(headers{1});
chr = fieldnames(x.win);

clear x

for k = 1:length(headers)
    y = eval(headers{k});
    Read = y.win.(chr{1})(2,:);
    FullGenome_Windows_top = y.win.(chr{1})(1,:);
    FullGenome_Windows_bottom = y.win.(chr{1})(1,:)+window-1;
    conden_len(1) = length(y.win.(chr{1}));
    for i = 2:length(chr)
        Read = [Read y.win.(chr{i})(2,:)];
        FullGenome_Windows_top = [FullGenome_Windows_top y.win.(chr{i})(1,:)];
        FullGenome_Windows_bottom = [FullGenome_Windows_bottom y.win.(chr{i})(1,:)+window-1];
        conden_len(i) = length(y.win.(chr{i}));
    end
    Reads(k,:) = Read./ReadSum(1,k);
    FullGenome_Windows = [FullGenome_Windows_top;FullGenome_Windows_bottom];
end
clear y k i FullGenome_Windows_bottom FullGenome_Windows_top Read

%% Set Parameters
									
win_thresh = 1000;									%The maximum distance betweeen 2 significant windows

reads_thresh = mean(mean(Reads(:,[IPidx controlidx])));   %Establish the reads threshold - Need to find a better statistical method to do this


%% Initialize cells
coordinates = cell(1,5);
sig_coordinates = cell(1,5);


%% Convert the data to the log scale

interg = normalize_fold_PreImmune(Reads,IPidx,controlidx,reads_thresh);

%% Determine threshold values via GUI and find significant/filtered windows

for k = 1:length(IPidx)
    thresh = 2*reads_thresh;
    
	%Filter against windows that have 0 enrichment or fall below the reads threshold
    i = interg(k,:) > 0 & (Reads(IPidx(k),:) + Reads(controlidx,:) > thresh);
    a = find(i);
    
    coordinates(k,:) = {Chr(i) FullGenome_Windows(:,i) Dist(i,:) ...
        Reads(IPidx(k),i) interg(k,i)};
    x.Chr = coordinates{k,1};
    x.Window = coordinates{k,2};
    x.Distance = coordinates{k,3};
    x.Reads = coordinates{k,4};
    x.Enrichment = coordinates{k,5};
    
    assignin('base',sprintf('%s_filtered',headers{IPidx(k)}),x);
    
    
    %Compute the zscore of the remaining Reads and identify significant
    %coordinates
    
    lam = poissfit(coordinates{k,5});
    Ab = 1 - poisscdf(coordinates{k,5},lam);
    sig = Ab < 1e-3;
    
    chrom = coordinates{k,1};
    windows = coordinates{k,2};
    distsubset = coordinates{k,3};
    Readsubset = coordinates{k,4};
    Enrichsubset = coordinates{k,5};
    sig_coordinates(k,:) = {chrom(sig) windows(1:2,sig) ...
        distsubset(sig,:) Readsubset(1,sig) Enrichsubset(1,sig)};
    y.Chr = sig_coordinates{k,1};
    y.Window = sig_coordinates{k,2};
    y.Distance = sig_coordinates{k,3};
    y.Reads = sig_coordinates{k,4};
    y.Enrichment = sig_coordinates{k,5};
    
    B = length(y.Distance);
    [Condensed Chr_conden TSS_conden] = neighbors_fromwindows_v3(B,win_thresh,y,knownGene,Reads(controlidx,a(sig)),chr);
    
    Thresh = (Condensed(:,7)*ReadSum(IPidx(k)) >= 15);
    
    %For each iteration save the structure to the appropriate IP
    assignin('base',sprintf('%s_significant',headers{IPidx(k)}),y);
    assignin('base',sprintf('%s_condensed',headers{IPidx(k)}),Condensed(Thresh,:));
    assignin('base',sprintf('%s_TSS',headers{IPidx(k)}),TSS_conden);
    assignin('base',sprintf('%s_Chr',headers{IPidx(k)}),Chr_conden(Thresh));
end 
clear chrsubset windows distsubset Readsubset control sig a i j k x y ...
    m sig_coordinates coordinates Logged B Enrichsubset Condensed ...
    win_thresh z ttry 