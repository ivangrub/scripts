dirs = {'Sample_TY012','Sample_TY013','Sample_TY014','Sample_TY015'};

for j = 1:length(dirs)
    cd(sprintf('/Volumes/Teppei/%s' ,dirs{j}));
    a = dir;
    [m,~] = size(a);
    for i = 1:m
        if strfind(a(i).name,'chip.') & strfind(a(i).name,'.bedgraph')
            fprintf(a(i).name)
            cd('~/Desktop')
            pf0_fromexpresswiggle(sprintf('/Volumes/Teppei/%s/%s',dirs{j},a(i).name),50,250,'~/Desktop/genome_mm9.mat')
        end	
    end
end