%Internal normalization

%% Initialize matrices

conden_len = zeros(1,35);

%% Process Structured Input
x = eval(headers{1});
chr = fieldnames(x.win);
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
    Reads(k,:) = Read;
    FullGenome_Windows = [FullGenome_Windows_top;FullGenome_Windows_bottom];
end
clear x y k i FullGenome_Windows_bottom FullGenome_Windows_top Read

%% Run Analysis

win_thresh = 1000;

sig_coordinates = cell(1,4);

for k = 1:length(headers)
    R = Reads(k,:);
    
    Win = FullGenome_Windows;
    Distance = Dist;
    
    Prob = 1 - poisscdf(R,mean(R));
    P = Prob <= 1e-32;
    
    sig_coordinates(k,:) = {Chr(P) Win(1:2,P) Distance(P,:) R(1,P)};
    y.Chr = sig_coordinates{k,1};
    y.Window = sig_coordinates{k,2};
    y.Distance = sig_coordinates{k,3};
    y.Reads = sig_coordinates{k,4};
    
    B = length(y.Distance);
    [Condensed Chr_conden TSS_conden] = neighbors_fromwindows_internal_SumReads(B,win_thresh,y,knownGene,chr);
    
    assignin('base',sprintf('%s_Poiss_Sig',headers{k}),y)
    assignin('base',sprintf('%s_Poiss_Conden',headers{k}),Condensed)
    assignin('base',sprintf('%s_TSS_Conden',headers{k}),TSS_conden);
    assignin('base',sprintf('%s_Chr_Conden',headers{k}),Chr_conden);
end