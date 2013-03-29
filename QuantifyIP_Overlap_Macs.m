%Calculate the Overlap of Macs peaks between MEF_TAF7 and MEF_TBP

clear,clc

IP = {'MEF_TBP' 'MEF_TAF7'};

name1 = 'TY_MEF_TBP_relativeto_MEF_PI_wTags_Macs.txt';
name2 = 'TY_MEF_TAF7_relativeto_MEF_PI_wTags_Macs.txt'; 

files = {name1 name2};
for i = 1:length(files)
    fid = fopen(files{i});
    data = textscan(fid,'%s%d%d%s%d%d%d%d%d%d','HeaderLines',1);
    fclose(fid);clear fid
    assignin('base',IP{i},data);
end


IND = 1;
for k = 1:length(MEF_TAF7{1,1})
    chr = strcmp(MEF_TBP{1,1},MEF_TAF7{1,1}(k));
    win = (MEF_TBP{1,2} >= MEF_TAF7{1,2}(k) & MEF_TBP{1,3} <= MEF_TAF7{1,3}(k)) | ...       %IP window fits inside of OS window
        (MEF_TBP{1,2} < MEF_TAF7{1,2}(k) & MEF_TBP{1,3} >= MEF_TAF7{1,2}(k)) | ...         %IP window straddles start of OS window
        (MEF_TBP{1,2} <= MEF_TAF7{1,3}(k) & MEF_TBP{1,3} > MEF_TAF7{1,3}(k)) | ...         %IP window staddles the end of the window
        (MEF_TBP{1,2} <= MEF_TAF7{1,2}(k) & MEF_TBP{1,3} >= MEF_TAF7{1,3}(k));              %OS window fits inside of the IP window
    Overlap = (chr & win);
    if sum(Overlap) == 1
        OSchr(IND) = MEF_TAF7{1,1}(k);
        OSstart(IND) = MEF_TAF7{1,2}(k);
        OSend(IND) = MEF_TAF7{1,3}(k);
        MEF_TAF7summit(IND) = MEF_TAF7{1,2}(k) + MEF_TAF7{1,5}(k);
        MEF_TBPchr(IND) = MEF_TBP{1,1}(chr & win);
        MEF_TBPstart(IND) = MEF_TBP{1,2}(chr & win);
        MEF_TBPend(IND) = MEF_TBP{1,3}(chr & win);
        MEF_TBPsummit(IND) = MEF_TBP{1,2}(chr & win) + MEF_TBP{1,5}(chr & win);
        OverlapTags(IND) = MEF_TBP{1,6}(chr & win);
        OverlapPvalue(IND) = MEF_TBP{1,7}(chr & win);
        OverlapEnrichment(IND) = MEF_TBP{1,8}(chr & win);
        IND = IND + 1;
        continue
    elseif sum(Overlap) > 1
        for RUN = 1:sum(Overlap)
            jj = find(Overlap);
            OSchr(IND) = MEF_TAF7{1,1}(k);
            OSstart(IND) = MEF_TAF7{1,2}(k);
            OSend(IND) = MEF_TAF7{1,3}(k);
            MEF_TAF7summit(IND) = MEF_TAF7{1,2}(k) + MEF_TAF7{1,5}(k);
            MEF_TBPchr(IND) = MEF_TBP{1,1}(jj(RUN));
            MEF_TBPstart(IND) = MEF_TBP{1,2}(jj(RUN));
            MEF_TBPend(IND) = MEF_TBP{1,3}(jj(RUN));
            MEF_TBPsummit(IND) = MEF_TBP{1,2}(jj(RUN)) + MEF_TBP{1,5}(jj(RUN));
            OverlapTags(IND) = MEF_TBP{1,6}(jj(RUN));
            OverlapPvalue(IND) = MEF_TBP{1,7}(jj(RUN));
            OverlapEnrichment(IND) = MEF_TBP{1,8}(jj(RUN));
            IND = IND + 1;
        end
        continue
    elseif sum(Overlap) == 0
        OSchr(IND) = MEF_TAF7{1,1}(k);
        OSstart(IND) = MEF_TAF7{1,2}(k);
        OSend(IND) = MEF_TAF7{1,3}(k);
        MEF_TAF7summit(IND) = MEF_TAF7{1,2}(k) + MEF_TAF7{1,5}(k);
        MEF_TBPchr(IND) = {'NaN'};
        MEF_TBPstart(IND) = 0;
        MEF_TBPend(IND) = 0;
        MEF_TBPsummit(IND) = 0;
        OverlapTags(IND) = 0;
        OverlapPvalue(IND) = 0;
        OverlapEnrichment(IND) = 0;
        IND = IND + 1;
        continue
    end
end   
fid = fopen(sprintf('OverlapPeaksFor_%s_relativeto_%s.txt',IP{1},IP{2}),'w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',sprintf('%s Chr',IP{2}),sprintf('%s Start',IP{2}),sprintf('%s End',IP{2}),...
    sprintf('%s Chr',IP{1}),sprintf('%s Start',IP{1}),sprintf('%s End',IP{1}),sprintf('%s Summit',IP{1}),sprintf('%s Tags',IP{1}), ...
    sprintf('%s Enrichment',IP{1}),sprintf('%s Pvalue',IP{1}));
for z = 1:length(OSchr)
    fprintf(fid,'%s\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n',OSchr{z},OSstart(z),OSend(z),MEF_TBPchr{z},MEF_TBPstart(z),MEF_TBPend(z),MEF_TBPsummit(z),...
        OverlapTags(z),OverlapPvalue(z),OverlapEnrichment(z));
end

NumPeaks = [length(MEF_TBP{1,1}) length(MEF_TAF7{1,1})];
cover = ~strcmp(MEF_TBPchr,'NaN');
Overlapped = sum(cover);
PercentagePeaks = sum(cover)/NumPeaks(2);

Distance = double(MEF_TAF7summit(cover)) - double(MEF_TBPsummit(cover));

[freq,xout] = hist(Distance,2000);
plot(xout,freq/length(Distance),'r')
xlabel('Distance (bp)')
ylabel('Fraction of Total Overlapping Peaks')
axis([-1000 1000 0 .5])

Summed = sum(freq(abs(xout)<=200))./length(Distance);

%saveas(gcf,sprintf('TeppeiFigures/SummitPlot_%s_relativeto_%s.eps',IP{1},IP{2}),'psc2');

%fclose(fid);clear fid IP IND OSchr OSstart OSend MEF_TBPchr MEF_TBPstart MEF_TBPend OverlapTags OverlapEnrichment OverlapPvalue Overlap ...
 %   MEF_TBPsummit RUN ans chr data files i jj k name1 name2 z

