%Take the base pair resolution structure and event the reads for the full coverage
%This will be specific for the 5kb region around the TSS

function avgreads = bp2TSS_wRNA(x,knownGene,read_length,win)

chr = fieldnames(x.bp);
bpres = zeros((win+read_length-1)*2+1,length(chr)+1);

for i = 1:length(knownGene)
	if cell2mat(knownGene(i,5)) == 1
		dist = x.bp.(knownGene{i,4}) - cell2mat(knownGene(i,6));
	else
		dist = (x.bp.(knownGene{i,4}) - cell2mat(knownGene(i,7)))*(-1);
	end
	f = find(abs(dist) < win);
	k = find(strcmp(knownGene(i,4),chr));
    k = k + 1;
    bpres(:,1) = (-win-read_length+1):(win+read_length-1);
    if ~isempty(f)
        for j = 1:length(f)
            bpres((dist(f(j)):dist(f(j))+read_length-2)+win,k) = bpres((dist(f(j)):dist(f(j))+read_length-2)+win,k) + 1;
        end
    end
end

avgreads = [bpres(:,1) mean(bpres(:,2:end),2)];

