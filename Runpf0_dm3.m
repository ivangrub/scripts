dirs = {'fly_pol2_300bp_neigh1_wB' 'fly_pol2_300bp_wB' 'fly_pol2_900bp_wB' 'fly_pol2_1000bp_neigh1_wB' 'fly_pol2_3000bp_wB' 'fly_pol2_3000bp_neigh1_wB' 'fly_pol2_9000bp_wB' 'fly_pol2_unique'};

for j = 1:length(dirs)
    cd(sprintf('~/Desktop/express_fly_pol2/%s' ,dirs{j}));
    a = dir;
    [m,~] = size(a);
    for i = 1:m
        if strfind(a(i).name,'chip.') & strfind(a(i).name,'.bedgraph')
            fprintf(a(i).name)
            cd('~/Desktop')
            pf0_fromexpresswiggle(sprintf('~/Desktop/express_fly_pol2/%s/%s',dirs{j},a(i).name),50,250,'~/Desktop/genome_dm3.mat')
        end	
    end
end