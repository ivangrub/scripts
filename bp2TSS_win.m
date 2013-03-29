%Look closer at a specific window with respect to the TSS

function gene_list = bp2TSS_win(x,knownGene,read_length,win)

chr = fieldnames(x.win);
bpres = zeros((win+read_length-1)*2+1,36);
genes = NaN(length(chr),length(knownGene));

for i = 1:length(knownGene)
	if cell2mat(knownGene(i,5)) == 1
		dist = x.bp.(knownGene{i,4}) - cell2mat(knownGene(i,6));
	else
		dist = (x.bp.(knownGene{i,4}) - cell2mat(knownGene(i,7)))*(-1);
	end
	f = find(abs(dist) >= win(1) & abs(dist) <= win(2));
	k = find(strcmp(knownGene(i,4),chr));
    k = k + 1;
    bpres(:,1) = (-win-read_length+1):(win+read_length-1);
    if ~isempty(f)
        for j = 1:length(f)
            bpres((dist(f(j)):dist(f(j))+read_length-2)+win,k) = bpres((dist(f(j)):dist(f(j))+read_length-2)+win,k) + 1;
            genes(j,i) = knownGene(i,2);
        end
    end
end

gene_list = sort(bpres,'descend');