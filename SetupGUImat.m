%Create structures of peaks, raw data for each IP to input into the GUI

clc

%Set directory where the .mat files for all of the IPs are
dir = '../ChIP-Seq Data/Teppei/';

person = 'TY_';
headers = {'TAF7' 'MEF_TAF7' 'MEF_PolII' 'MEF_TBP' 'MEF_PI' 'PI'};
peaks = {'Poisson_Conden' 'Lin_Conden' 'Macs'};

for i = 1:length(headers)
    x = load(sprintf('%s%s%s_SumReads.mm9.mat',dir,person,headers{i}));
    assignin('base',sprintf('%s',headers{i}),x.(headers{i}).win);
    clear x
end
clear dir i

control = {'PI' 'MEF_PI' 'MEF_PI' 'MEF_PI'};
for i = 1:4
    for j = 1:length(peaks)
        if ~strcmp(peaks{j},'Poisson_Conden')
            y = eval(headers{i});
            x = importdata(sprintf('%s%s_relativeto_%s_%s.txt',person,headers{i},control{i},peaks{j}));
            y.(peaks{j}) = x.data;
            y.(sprintf('%s_text',peaks{j})) = x.textdata;
            assignin('base',sprintf('%s',headers{i}),y);
            clear x y
        else
            y = eval(headers{i});
            x = importdata(sprintf('%s%s_%s.txt',person,headers{i},peaks{j}));
            y.(peaks{j}) = x.data;
            y.(sprintf('%s_text',peaks{j})) = x.textdata;
            assignin('base',sprintf('%s',headers{i}),y);
            clear x y
        end
    end
end

clear headers i j peaks person control

load knownGene.mm9.mat

clear CpG