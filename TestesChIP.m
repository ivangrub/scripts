
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
		coordind = ceil(ind*50);
		chromind = peak(j,1) == known(:,1);
		chrknown = known(chromind,:)
		midpoint = mean(peak(2),peak(3));
		promoter = midpoint >= chrknown(:,2)-500 | midpoint <= chrknown(:,2)+500;
		proximal = midpoint >= chrknown(:,2)-10000 | midpoint <= chrknown(:,2)+10000;
		distal = midpoint >= chrknown(:,2)-50000 | midpoint <= chrknown(:,2)+50000;
		if sum(promoter) > 0
			region = 'promoter';
		else if sum(proximal) > 0
			region = 'proximal';
		else if sum(proximal) > 0
			region = 'distal';
		else
			region = 'intergrenic';
		end

		[tssind,closesttss] = min(abs(midpoint - chrknown(:,2)));
		tssgenename = chrknown(tssind,1);

		genelist = cell(1,3);
		genedist = zeros(1,3);
		[index,distance] = sort(abs(midpoint - chrknown(:,2)));
		for closgene = 1:3
			genelist{closgene} = chrknown(index(closgene),1);
			genedist(closgene) = distance(closgene);
			tommytssdistance = midpoint - chrknown(index(closgene),2);
			tommytesdistance  = midpoint - chrknown(index(closgene),3);
			if tommytssdistance < 0
				genedist(closgene) = genedist(closgene)*-1;
			else if tommytssdistance > 0 && tommytesdistance < 0
				genedist(closgene) = 0;
			else
				continue
			end
		end

		range = ceil(peak(2)/50):ceil(peak(3)/50);
		chrom = peak(j,1) == peak(:,1);
		[ind,val] = max(ip(chrom,range));

		for k = 1:length(IP)
			if ~strcmp(IP{i},IP{k})
				y  = load(sprintf('%s.mat',IP{k}));
			else 
				continue
			end
			maxenrich = max(y(chrom,range))
		end
	end

end