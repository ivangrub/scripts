clear

fid = fopen('transcripts.gtf','r');

A = textscan(fid,'%s%s%s%d%d%d%s%s%s%s%s%s%s%s%s%s','delimiter','\t');
fclose(fid);clear fid

TranscriptOnly = strcmp(A{1,3},'transcript');

FPKM = A{1,11}(TranscriptOnly);
FPKM = strrep(FPKM,'FPKM "','');
FPKM = strrep(FPKM,'"','');
FPKM = str2double(FPKM);

GI = A{1,10}(TranscriptOnly);
GI = strrep(GI,'transcript_id "','');
GI = strrep(GI,'"','');

discFPKM = FPKM(FPKM < 1);
discGI = GI(FPKM < 1);

fid = fopen('NotTranscribed.ES.txt','w');
for i = 1:length(discFPKM)
     fprintf(fid,'%s\t%d\n',char(discGI{i}),discFPKM(i));
end
fclose(fid);clear fid

newFPKM = FPKM(FPKM >= 1);
GI = GI(FPKM >= 1);
[sortedFPKM,j] = sort(newFPKM,'ascend');
sortedGI = GI(j);

bin = round(linspace(1,length(sortedFPKM),5));
bins = [100 75 50 25 0];
for i = 1:length(bins)-1
    fid = fopen(sprintf('Between_%s_and_%s_Percent.RNASeqTranscripts.ES.txt',num2str(bins(i)),num2str(bins(i+1))),'w');
    
    for j = bin(i):bin(i+1)
        fprintf(fid,'%s\t%d\n',char(sortedGI{j}),sortedFPKM(j));
    end
    fclose(fid);clear fid
end
    