
IP = {'TAF7L' 'TBP' 'PolII' 'AR' 'TAF7' 'PI'};
IPlist = IP;

geneid = knownGene(:,1);
for i = 1:length(knownGene)
	if strcmp(knownGene{i,4},'+')
		continue
	else
		st = knownGene{i,3};
        endi = knownGene{i,2};
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

fid = fopen('Testes_geneChIP.txt','w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',(sprintf('%s_promoter',IP{1}),sprintf('%s_proximal',IP{1}),...
	sprintf('%s_distal',IP{1}),sprintf('%s_promoter',IP{2}),sprintf('%s_proximal',IP{2}),sprintf('%s_distal',IP{2}),...
	sprintf('%s_promoter',IP{3}),sprintf('%s_proximal',IP{3}),sprintf('%s_distal',IP{3}),...
	sprintf('%s_promoter',IP{4}),sprintf('%s_proximal',IP{4}),sprintf('%s_distal',IP{4}),...
	sprintf('%s_promoter',IP{5}),sprintf('%s_proximal',IP{5}),sprintf('%s_distal',IP{5}),...
	sprintf('%s_promoter',IP{6}),sprintf('%s_proximal',IP{6}),sprintf('%s_distal',IP{6}))
for i = 1:length(enrich)
	fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',(enrich(i,1),enrich(i,2),enrich(i,3),...
		enrich(i,4),enrich(i,5),enrich(i,6),enrich(i,7),enrich(i,8),enrich(i,9),enrich(i,9),enrich(i,10),enrich(i,11),enrich(i,12),...
		enrich(i,13),enrich(i,14),enrich(i,15),enrich(i,16),enrich(i,17));
end
fclose(fid),clear fid
for i = 1:length(IP)
	ip = eval(IP{i}); %load(sprintf('%s.mat',IP{i}));
	fid = fopen(sprintf('%s.filtered.peaks.bed',IP{i}));
	peak = textscan(fid,'%s%d%d','delimiter','\t');
	fclose(fid);clear fid;
	otherIP = find(IP(~IP{i}));
	fid = fopen(sprintf('%s_peakcentric.txt',IP{i});
	fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',('Chromosome','Left','Right','Cis-Reg',...
		'Closest Gene','Distance','2nd Closest','Distance','3rd Closest','Distance',IP{otherIP(1)},...
		IP{otherIP(2)},IP{otherIP(3)},IP{otherIP(4)},IP{otherIP(5})
		
	maxenrich = zeros(1,length(IP)-1)
	genesect = zeros(length(peak),length(IP)+4);
	for j = 1:length(peak)
		coordind = ceil(ind*50);
		chromind = peak(j,1) == known(:,1);
		chrknown = known(chromind,:)
		gidknown = geneid(chromind,:)

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
			genelist{closgene} = gidknown(index(closgene),1);
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

		for k = 1:length(IPlist)
			if ~strcmp(IP{i},IPlist{k})
				y  = eval(IP{i}); %load(sprintf('%s.mat',IPlist{k}));
			else 
				continue
			end
			maxenrich(k) = max(y(chrom,range))
		end
		
		fprintf(fid,'%s\t%d\t%d\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n',(peak(j,1),peak(j,2),peak(j,3),region,genelist(1),...
			genedist(1),genelist(2),genedist(2),genelist(3),genedist(3),maxenrich(1),maxenrich(2),maxenrich(3),maxenrich(4),maxenrich(5)));
		
	end
	fclose(fid);clear fid
end