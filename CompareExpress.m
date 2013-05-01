conditions = {'pol2_100bp_neigh1_wB','pol2_100bp_neigh1_wB_samp','pol2_100bp_samp_wB','pol2_300bp_neigh1_wB','pol2_300bp_wB_samp','pol2_300bp_wB','noexpress'};
names = cell(1,length(conditions));
for z = 1:length(conditions)

    cd(sprintf('/Volumes/Genomic1/IG_express/fly_pol2_k100_fakerepeat/%s',conditions{z}));
    a = dir;

    bundle = 1000;
    [m,~] = size(a);
    
    for i = 1:m
        if ~isempty(strfind(a(i).name,'chip.') & strfind(a(i).name,'.bedgraph'))
            fprintf('On %s\n',a(i).name)
            fid = fopen(a(i).name,'r');
            if z == 1
                [s,w] = system(sprintf('wc -l %s',a(i).name));
                [len,~] = strread(w,'%d%s');
                FullCount = zeros(length(conditions),len);
                A = zeros(1,len);
            end
            k = 1;
            while ~feof(fid)
                B = textscan(fid,'%s%d%f',bundle,'HeaderLines',1);
                if k + bundle - 1 > len
                    A(1,k:k + length(B{3})-1) = B{3};
                else
                    A(1,k:k+bundle-1) = B{3};
                    k = k + bundle;
                end   
            end
            fclose(fid);clear fid B
        end
    end
    FullCount(z,:) = A;
    names(z) = conditions(z);
end