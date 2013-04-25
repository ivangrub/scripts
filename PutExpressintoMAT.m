
conditions = {'pol2_100bp_neigh1_wB','pol2_300bp_wB','pol2_300bp_neigh1_wB','pol2_900bp_wB','pol2_9000bp_wB'};
thisfile = 'ReadProbability.txt';

for j = 1:length(conditions)
    cd(sprintf('/Volumes/Genomic1/IG_express/%s' ,conditions{j}));
    a = dir;

    bundle = 1000;
    [m,~] = size(a);
    j = 1;
    names = cell(1,m);
    for i = 1:m
        if ~isempty(strfind(a(i).name,thisfile))
            fprintf('On %s\n',a(i).name)
            fid = fopen(a(i).name,'r');
            if j == 1
                [s,w] = system(sprintf('wc -l %s',a(i).name));
                [len,~] = strread(w,'%d%s');
                A = zeros(m,len);
            end
            k = 1;
            while ~feof(fid)
                B = textscan(fid,'%s%d%f',bundle,'HeaderLines',1);
                if k + bundle - 1 > len
                    A(j,k:k + length(B{3})-1) = B{3};
                else
                    A(j,k:k+bundle-1) = B{3};
                    k = k + bundle;
                end   
            end
            names{j} = a(i).name;
            j = j + 1;
            fclose(fid);clear fid B
        end
    end
    FullCount = A(1:j-1,:);
    names = names(1:j-1);
end