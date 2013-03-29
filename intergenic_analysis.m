%Find the Distance to nearest TSS

[m,~] = size(Chr);
Dist = zeros(m,2);
chr_num = strrep(chr_num,'chr','');

Win = mean(FullGenome_Windows,2);

%Find the minimum distance to the closest TSS for each of the enriched
%intergenic regions with respect to each IP

k = 1;
for j = 1:length(chr_num)
    chr = strcmp(Chr,chr_num(j));
    chr_TSS = strcmp(chr_num(j),knownGene(:,4));
    o = find(chr_TSS);
    known2 = cell2mat(knownGene(o,6));
    for i = 1:sum(chr)
        Dist2Win = Win(k,1) - known2;
        [dist2tss,b] = min(abs(Dist2Win));
%         if dist2tss ~= diff(b)
%             Diff = diff(b);
%         end
        Dist(k,1:2) = [Win(k,1) dist2tss];
        k = k + 1;
    end
end
