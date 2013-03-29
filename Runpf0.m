dirs = {'sim_300bp_neigh1_wB' 'sim_300bp_wB' 'sim_900bp_wB' 'sim_1000bp_neigh1_wB' 'sim_3000bp_wB' 'sim_3000bp_neigh1_wB' 'sim_9000bp_wB' 'sim_unique'};

for j = 1:length(dirs)
    cd(sprintf('~/Desktop/express_sim/%s' ,dirs{j}));
    a = dir;
    [m,~] = size(a);
    for i = 1:m
        if strfind(a(i).name,'chip.') & strfind(a(i).name,'.bedgraph')
            fprintf(a(i).name)
            cd('~/Desktop')
            pf0_fromexpresswiggle(sprintf('~/Desktop/express_sim/%s/%s',dirs{j},a(i).name),50,250,'~/Desktop/genome_mm9.mat')
        end	
    end
end