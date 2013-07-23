type = {'CSEM' 'Express' 'NoExpress' 'Unique'};
IP = {'Pol2' 'TBP' 'PI'};
annot = {'refgene' 'ensgene' 'knownGene' 'cufflinks_denovo'};


for anno = 1:length(annot)
    k = 1;
    disp(sprintf('Gene centric using %s annotation',annot{anno}))
    x = load(sprintf('%s.mat',annot{anno}));
    knownGene = x.(annot{anno});
    geneid = knownGene(:,4);
    genechrom = knownGene(:,1);
    known = zeros(length(knownGene),3);
    for i = 1:length(knownGene)
        if strcmp(knownGene{i,6},'+')
            st = knownGene{i,2};
            endi = knownGene{i,3};
        else
            st = knownGene{i,2};
            endi = knownGene{i,3};
        end
        known(i,2) = st;
        known(i,3) = endi;
    end
    
    enrich = zeros(length(known),3*length(IP)*length(type));
    
    for typ = 1:length(type)

        for i = 1:length(IP)
            y = load(sprintf('%s_%s.mat',type{typ},IP{i}));
            for j = 1:length(known)
                if strfind(genechrom{j},'random') > 0
                    continue
                end
                
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
    end
    
    fid = fopen(sprintf('GeneChIP_%s.txt',annot{anno}),'w');
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Gene ID',sprintf('%s_%s_promoter',type{1},IP{1}),sprintf('%s_%s_proximal',type{1},IP{1}),...
        sprintf('%s_%s_distal',type{1},IP{1}),sprintf('%s_%s_promoter',type{1},IP{2}),sprintf('%s_%s_proximal',type{1},IP{2}),sprintf('%s_%s_distal',type{1},IP{2}),...
        sprintf('%s_%s_promoter',type{1},IP{3}),sprintf('%s_%s_proximal',type{1},IP{3}),sprintf('%s_%s_distal',type{1},IP{3}),...
        sprintf('%s_%s_promoter',type{2},IP{1}),sprintf('%s_%s_proximal',type{2},IP{1}),sprintf('%s_%s_distal',type{2},IP{1}),...
        sprintf('%s_%s_promoter',type{2},IP{2}),sprintf('%s_%s_proximal',type{2},IP{2}),sprintf('%s_%s_distal',type{2},IP{2}),...
        sprintf('%s_%s_promoter',type{2},IP{3}),sprintf('%s_%s_proximal',type{2},IP{3}),sprintf('%s_%s_distal',type{2},IP{3}),...
        sprintf('%s_%s_promoter',type{3},IP{1}),sprintf('%s_%s_proximal',type{3},IP{1}),sprintf('%s_%s_distal',type{3},IP{1}),...
        sprintf('%s_%s_promoter',type{3},IP{2}),sprintf('%s_%s_proximal',type{3},IP{2}),sprintf('%s_%s_distal',type{3},IP{2}),...
        sprintf('%s_%s_promoter',type{3},IP{3}),sprintf('%s_%s_proximal',type{3},IP{3}),sprintf('%s_%s_distal',type{3},IP{3}),...
        sprintf('%s_%s_promoter',type{4},IP{1}),sprintf('%s_%s_proximal',type{4},IP{1}),sprintf('%s_%s_distal',type{4},IP{1}),...
        sprintf('%s_%s_promoter',type{4},IP{2}),sprintf('%s_%s_proximal',type{4},IP{2}),sprintf('%s_%s_distal',type{4},IP{2}),...
        sprintf('%s_%s_promoter',type{4},IP{3}),sprintf('%s_%s_proximal',type{4},IP{3}),sprintf('%s_%s_distal',type{4},IP{3}));
    
    for i = 1:length(enrich)
        if isnumeric(geneid{i})
            id = num2str(geneid{i});
        else
            id = geneid{i};
        end
        fprintf(fid,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',char(id),enrich(i,1),enrich(i,2),enrich(i,3),...
            enrich(i,4),enrich(i,5),enrich(i,6),enrich(i,7),enrich(i,8),enrich(i,9),enrich(i,10),enrich(i,11),enrich(i,12),...
            enrich(i,13),enrich(i,14),enrich(i,15),enrich(i,16),enrich(i,17),enrich(i,18),enrich(i,19),enrich(i,20), ...
            enrich(i,21),enrich(i,22),enrich(i,23),enrich(i,24),enrich(i,25),enrich(i,26),enrich(i,27),enrich(i,28),enrich(i,29),...
            enrich(i,30),enrich(i,31),enrich(i,32),enrich(i,33),enrich(i,34),enrich(i,35),enrich(i,36));
    end
    fclose(fid);clear fid
end

for i = 1:length(type)
    for j = 1:length(IP)
        x = load(sprintf('%s_%s.mat',type{i},IP{j}));
        assignin('base',sprintf('%s_%s',type{i},IP{j}),x);
    end
end

for anno = 1:length(annot)
    disp(sprintf('Peak centric using %s annotation',annot{anno}))
    x = load(sprintf('%s.mat',annot{anno}));
    knownGene = x.(annot{anno});
    geneid = knownGene(:,4);
    genechrom = knownGene(:,1);
    known = zeros(length(knownGene),4);
    for i = 1:length(knownGene)
        if strcmp(knownGene{i,6},'+')
            st = knownGene{i,2};
            endi = knownGene{i,3};
            str = 1;
        else
            st = knownGene{i,2};
            endi = knownGene{i,3};
            str = -1;
        end
        known(i,2) = st;
        known(i,3) = endi;
        known(i,4) = str;
    end
    for typ = 1:length(type)
        for i = 1:length(IP)-1
            disp(sprintf('Working on %s_%s',type{typ},IP{i}))
            fid = fopen(sprintf('Peaks/Filtered/%s_%s.filtered.bed',type{typ},IP{i}));
            peaks = textscan(fid,'%s%d%d%f','delimiter','\t');
            fclose(fid);clear fid;
            
            fid = fopen(sprintf('%s_%s_peakcentric.txt',type{typ},IP{i}),'w');
            fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','Chromosome','Left','Right','Cis-Reg','Closest TSS','Distance',...
                'Closest Gene','Distance','2nd Closest','Distance','3rd Closest','Distance','Express_Pol2',...
                'Express_TBP','Express_PI','NoExpress_Pol2','NoExpress_TBP','NoExpress_PI','Unique_Pol2','Unique_TBP','Unique_PI',...
                'CSEM_Pol2','CSEM_TBP','CSEM_PI');
            
            maxenrich = zeros(1,12);
            genesect = zeros(length(peaks),length(type)*length(IP)+4);
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
                %   region = 'promoter';
                %elseif sum(proximal) > 0
                %   region = 'proximal';
                %elseif sum(proximal) > 0
                %   region = 'distal';
                %else
                %   region = 'intergrenic';
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
                
                
                maxenrich(1) = max(Express_Pol2.('chip').(chrom{1})(range));
                maxenrich(2) = max(Express_TBP.('chip').(chrom{1})(range));
                maxenrich(3) = max(Express_PI.('chip').(chrom{1})(range));
                maxenrich(4) = max(NoExpress_Pol2.('chip').(chrom{1})(range));
                maxenrich(5) = max(NoExpress_TBP.('chip').(chrom{1})(range));
                maxenrich(6) = max(NoExpress_PI.('chip').(chrom{1})(range));
                maxenrich(7) = max(Unique_Pol2.('chip').(chrom{1})(range));
                maxenrich(8) = max(Unique_TBP.('chip').(chrom{1})(range));
                maxenrich(9) = max(Unique_PI.('chip').(chrom{1})(range));
                maxenrich(10) = max(CSEM_Pol2.('chip').(chrom{1})(range));
                maxenrich(11) = max(CSEM_TBP.('chip').(chrom{1})(range));
                maxenrich(12) = max(CSEM_PI.('chip').(chrom{1})(range));
                
                if iscellstr(genelist{1}) == 0
                    for list = 1:length(genelist)
                        genelist{list} = num2str(cell2mat(genelist{list}));
                    end
                end
                
                fprintf(fid,'%s\t%d\t%d\t%s\t%s\t%s\t%s\t%d\t%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',char(peaks{1,1}(j)),peak(j,1),peak(j,2),region,tssgenename,num2str(closesttss),char(genelist{1}),...
                    genedist(1),char(genelist{2}),genedist(2),char(genelist{3}),genedist(3),maxenrich(1),maxenrich(2),maxenrich(3),maxenrich(4),maxenrich(5),maxenrich(6),...
                    maxenrich(7),maxenrich(8),maxenrich(9),maxenrich(10),maxenrich(11),maxenrich(12));
                
            end
            fclose(fid);clear fid
        end
    end
end