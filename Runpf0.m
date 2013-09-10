dirs = {'Desktop'};

for j = 1:length(dirs)
    cd(sprintf('~/%s' ,dirs{j}));
    a = dir;
    [m,~] = size(a);
    for i = 1:m
        if strfind(a(i).name,'chip.ES_TBP') & strfind(a(i).name,'.bedgraph')
            fprintf(a(i).name)
            if ~strcmp(a(i).name,'chip.CSEM_ES_Pol2.k100.bedgraph')
                continue
            end
            cd('~/Desktop')
            pf0_fromexpresswiggle(sprintf('~/%s/%s',dirs{j},a(i).name),50,250,'~/Desktop/genome_mm9.mat')
        end	
    end
end