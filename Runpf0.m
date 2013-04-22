dirs = {'ES_PI','ES_Pol2','ES_PI_150bp_neigh1_wB','ES_Pol2','ES_Pol2_150bp_neigh1_wB','ES_Pol2_neigh1_wB_samp','ES_TBP','ES_TBP_150bp_neigh1_wB'};

for j = 1:length(dirs)
    cd(sprintf('/Volumes/Genomic1/IG_express/%s' ,dirs{j}));
    a = dir;
    [m,~] = size(a);
    for i = 1:m
        if strfind(a(i).name,'chip.') & strfind(a(i).name,'.bedgraph')
            fprintf(a(i).name)
            cd('~/Desktop')
            pf0_fromexpresswiggle(sprintf('/Volumes/Genomic1/IG_express/%s/%s',dirs{j},a(i).name),50,250,'~/Desktop/genome_mm9.mat')
        end	
    end
end