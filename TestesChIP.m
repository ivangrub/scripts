
IP = {'TAF7L' 'TBP' 'PolII' 'AR' 'TAF7' 'PI'};

fid = fopen('knownGene.txt');

known = textscan(fid,'%s%d%d%s','delimiter','\t');
fclose(fid);clear fid
for i = 1:length(known)
	if strcmp(known{i,4},'+')
		continue
	else
		st = known{i,3};
        endi = known{i,2};
	end
	known(i,2) = st;
	known(i,3) = endi;
end

enrich = zeros(length(known),3*length(IP));

k = 1;
for i = 1:length(IP)
	y = load(sprintf('%s.mat',IP{i}));
	for j = 1:length(known)
		st = known{j,2};
		coord = ceil(st/50);
		ind = y(1) == known{j,1};
		prom = 500/50;
		prox = 10000/50;
		dist = 50000/50;
		enrich(j,k) = max(y(ind,2)-prom:y(ind,2)+prom);
		enrich(j,k+1) = max(max(y(ind,2)-prox:y(ind,2)-prom),max(y(ind,2)+prom:y(ind,2)+prox));
		enrich(j,k+2) = max(max(y(ind,2)-dist:y(ind,2)-prox),max(y(ind,2)+prox:y(ind,2)+dist));
	end 
	k = k + 3;
end

for i = 1:length(IP)
	ip = load(sprintf('%s.mat',IP{i}));
	fid = fopen(sprintf('%s.filtered.peaks.bed',IP{i}));
	peak = textscan(fid,'%s%d%d','delimiter','\t');
	fclose(fid);clear fid;
	genesect = zeros(length(peak),length(IP)+4);
	for j = 1:length(peak)
		range = ceil(peak(2)/50):ceil(peak(3)/50);
		chrom = peak(j,1) == peak(:,1);
		[ind,val] = max(ip(chrom,range));
		coordind = ceil(ind*50);
		chromind = peak(j,1) == known(:,1);
		[index,diff] = min(abs(ceil(known(chromind,2)/50) - coordind));
		
	end
	for k = 1:length(IP)
		if ~strcmp(IP{i},IP{k})
			y  = load(sprintf('%s.mat',IP{k}));
		else 
			continue
		end

	end
end