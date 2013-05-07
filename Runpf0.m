dirs = {'ES_TAF1_unique' 'ES_TAF1_multi' 'ES_TAF7_unique' 'ES_TAF7_multi' 'MEF_TAF7_unique' 'MEF_TAF7_multi'};

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