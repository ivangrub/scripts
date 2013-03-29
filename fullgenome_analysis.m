%% Full Genome Data Normalization and Analysis

clc

%% Set Parameters

win = (FullGenome_Windows(1,2)-FullGenome_Windows(1,1))*25-24;		%Calculate the size of the window 
reads_thresh = 2;													%Establish the reads threshold - Need to find a better statistical method to do this
win_thresh = 1000;													%The maximum distance betweeen 2 significant windows
control = find(strcmp('PI',headers)); 

if control > length(headers)
    error('Where is the PI control in the headers?')
end

%% Initialize cells
coordinates = cell(1,4);
sig_coordinates = cell(1,4);


%% Convert the data to the log scale
%Logged = log2(Reads + 1);
interg = normalize_fold(Reads,control,reads_thresh,'Lin');

%% Determine threshold values via GUI and find significant/filtered windows

for k = 1:(control-1)
    thresh = 10;
    
	%Filter against windows that have 0 enrichment or fall below the reads threshold
    i = interg(:,k) ~= 0 & (Reads(:,k) + Reads(:,control) >= thresh);
    a = find(i);
    coordinates(k,:) = {[FullGenome_Windows(a,1)*25-24 FullGenome_Windows(a,2)*25] ...
        Dist(a,:) Reads(a,k) interg(a,k)};
    x.Window = coordinates{k,1};
    x.Distance = coordinates{k,2};
    x.Reads = coordinates{k,3};
    x.Enrichment = coordinates{k,4};
    
    assignin('base',sprintf('%s_filtered',headers{k}),x);
    
    
    %Compute the zscore of the remaining Reads and identify significant
    %coordinates
    Ab = zscore(log2(coordinates{k,4}));
    sig = find(Ab >= 3);
    
    windows = coordinates{k,1};
    distsubset = coordinates{k,2};
    Readsubset = coordinates{k,3};
    Enrichsubset = coordinates{k,4};
    sig_coordinates(k,:) = {windows(sig,1:2) ...
        distsubset(sig,:) Readsubset(sig,1) Enrichsubset(sig,1)};
    y.Window = sig_coordinates{k,1};
    y.Distance = sig_coordinates{k,2};
    y.Reads = sig_coordinates{k,3};
    y.Enrichment = sig_coordinates{k,4};
    
     B = length(y.Distance);
     Condensed = neighbors_fromwindows_v2(B,win_thresh,win,y,knownGene,Reads(a(sig),control));
     
    %For each iteration save the structure to the appropriate IP
    assignin('base',sprintf('%s_significant',headers{k}),y);
    assignin('base',sprintf('%s_condensed',headers{k}),Condensed); 
end 
clear chrsubset windows distsubset Readsubset control Ab sig a i j k x y ...
    m sig_coordinates coordinates interg Logged B Enrichsubset Condensed ...
    win_thresh z ttry