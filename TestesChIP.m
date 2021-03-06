
IP = {'TAF7L' 'TBP' 'PolII' 'TAF7' 'AR'};
IPlist = IP;

geneid = knownGene(:,1);
genechrom = knownGene(:,2);
known = zeros(length(knownGene),4);
for i = 1:length(knownGene)
    if strcmp(knownGene{i,3},'+')
        st = knownGene{i,4};
        endi = knownGene{i,5};
        str = 1;
    else
        st = knownGene{i,5};
        endi = knownGene{i,4};
        str = -1;
    end
	known(i,2) = st;
	known(i,3) = endi;
	known(i,4) = str;
end

enrich = zeros(length(known),3*length(IP)-1);

k = 1;
for i = 1:length(IP)
	y = load(sprintf('%s.mat',IP{i}));
	for j = 1:length(known)
		edge = zeros(1,3);

		st = known(j,2);
		coord = ceil(st/50);
        
		ind = y.('chip').(genechrom{j});
		prom = 500/50;
		prox = 10000/50;
		dist = 50000/50;
        
        promreg = [coord-prom:coord+prom];
        XXprom = promreg > 0 & promreg < length(ind);
        if sum(XXprom) > 0
            enrich(j,k) = max(ind(promreg(XXprom)));
        else
            enrich(j,k) = NaN;
        end
        
        proxreg = [coord-prox:coord-prom,coord+prom:coord+prox];
        XXprox = proxreg > 0 & proxreg < length(ind);
        if sum(XXprox) > 0
            enrich(j,k+1) = max(ind(proxreg(XXprox)));
        else
            enrich(j,k+1) = NaN;
        end
        
        distreg = [coord-dist:coord-prox,coord+prox:coord+dist];
        XXdist = distreg > 0 & distreg < length(ind);
        if sum(XXdist) > 0
            enrich(j,k+2) = max(ind(distreg(XXdist)));
        else
            enrich(j,k+2) = NaN;
        end
        
	end 
	k = k + 3;
end

fid = fopen('Testes_geneChIP.txt','w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Gene ID','Gene Description',sprintf('%s_promoter',IP{1}),sprintf('%s_proximal',IP{1}),...
	sprintf('%s_distal',IP{1}),sprintf('%s_promoter',IP{2}),sprintf('%s_proximal',IP{2}),sprintf('%s_distal',IP{2}),...
	sprintf('%s_promoter',IP{3}),sprintf('%s_proximal',IP{3}),sprintf('%s_distal',IP{3}),...
	sprintf('%s_promoter',IP{4}),sprintf('%s_proximal',IP{4}),sprintf('%s_distal',IP{4}),...
	sprintf('%s_promoter',IP{5}),sprintf('%s_proximal',IP{5}),sprintf('%s_distal',IP{5}));
for i = 1:length(enrich)
    index = 0;
    if sum(index) == 0
       genename = {'unknown'}; 
    else
       genename = gidadesc(index,2);
    end
	fprintf(fid,'%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',char(geneid{i}),char(genename),enrich(i,1),enrich(i,2),enrich(i,3),...
		enrich(i,4),enrich(i,5),enrich(i,6),enrich(i,7),enrich(i,8),enrich(i,9),enrich(i,10),enrich(i,11),enrich(i,12),...
        enrich(i,13),enrich(i,14),enrich(i,15));
end
fclose(fid);clear fid
TAF7L = load('TAF7L.mat');
TBP = load('TBP.mat');
Pol2 = load('PolII.mat');
TAF7 = load('TAF7.mat');
AR = load('AR.mat');
PI = load('IgG.mat');

for i = 1:length(IP)
   
	fid = fopen(sprintf('%s.filtered.peaks.bed',IP{i}));
	peaks = textscan(fid,'%s%d%d%f','delimiter','\t');
	fclose(fid);clear fid;
	
	fid = fopen(sprintf('%s_peakcentric.txt',IP{i}),'w');
	fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Chromosome','Left','Right','Cis-Reg','Closest TSS','Distance',...
		'Closest Gene','Distance','2nd Closest','Distance','3rd Closest','Distance',char(IP{1}),...
		char(IP{2}),char(IP{3}),char(IP{4}),char(IP{5}),'IgG');
		
	maxenrich = zeros(1,length(IP)-1);
	genesect = zeros(length(peaks),length(IP)+4);
    peak = [peaks{2},peaks{3}];
	for j = 1:length(peaks{1})
		
		chromind = strcmp(peaks{1,1}(j),genechrom(:,1));
		chrknown = known(chromind,:);
		gidknown = geneid(chromind,:);

		midpoint = mean(peak(j,1),peak(j,2));
		promoter = midpoint >= chrknown(:,2)-500 | midpoint <= chrknown(:,2)+500;
		proximal = midpoint >= chrknown(:,2)-10000 | midpoint <= chrknown(:,2)+10000;
		distal = midpoint >= chrknown(:,2)-50000 | midpoint <= chrknown(:,2)+50000;

		%if sum(promoter) > 0
		%	region = 'promoter';
        %elseif sum(proximal) > 0
		%	region = 'proximal';
        %elseif sum(proximal) > 0
		%	region = 'distal';
		%else
		%	region = 'intergrenic';
		%end

		[closesttss,tssind] = min(abs(midpoint - chrknown(:,2)));
		tssgenename = gidknown{tssind};

		genelist = cell(1,3);
		genedist = zeros(1,3);
		[distance,index] = sort(abs(midpoint - chrknown(:,2)));
        
        if distance(1) <= 500
			region = 'promoter';
		elseif distance(1) > 500 && distance(1) <= 5000
			region = 'proximal';
		elseif distance(1) > 5000 && distance(1) <= 25000
			region = 'distal';
		else
			region = 'intergenic';
        end
        
		for closgene = 1:3
			genelist{closgene} = gidknown(index(closgene),1);
			genedist(closgene) = distance(closgene);
			
			tommytssdistance = midpoint - chrknown(index(closgene),2);
			tommytesdistance  = midpoint - chrknown(index(closgene),3);
			if tommytssdistance < 0
				genedist(closgene) = genedist(closgene)*-1*chrknown(index(closgene),4);
            elseif tommytssdistance > 0 && tommytesdistance < 0
				genedist(closgene) = 0;
			elseif tommytssdistance > 0 && tommytesdistance > 0
				genedist(closgene) = tommytesdistance*chrknown(index(closgene),4);
			end
		end

		range = ceil(peak(j,1)/50):ceil(peak(j,2)/50);
		chrom = peaks{1,1}(j);

		
        maxenrich(1) = max(TAF7L.('chip').(chrom{1})(range));
        maxenrich(2) = max(TBP.('chip').(chrom{1})(range));
        maxenrich(3) = max(Pol2.('chip').(chrom{1})(range));
        maxenrich(4) = max(TAF7.('chip').(chrom{1})(range));
        maxenrich(5) = max(AR.('chip').(chrom{1})(range));
        maxenrich(6) = max(PI.('chip').(chrom{1})(range));
		fprintf(fid,'%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',chrom{1},peak(j,1),peak(j,2),region,tssgenename,num2str(closesttss),char(genelist{1}),...
			genedist(1),char(genelist{2}),genedist(2),char(genelist{3}),genedist(3),maxenrich(1),maxenrich(2),maxenrich(3),maxenrich(4),maxenrich(5),maxenrich(6));
		
	end
	fclose(fid);clear fid
end