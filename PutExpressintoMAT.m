
conditions = {'pol2_100bp_neigh1_wB','pol2_300bp_wB','pol2_300bp_neigh1_wB','pol2_900bp_wB','pol2_9000bp_wB'};
thisfile = 'ReadProbability.txt';

for z = 1:length(conditions)
    cd(sprintf('/Volumes/Genomic1/IG_express/fly_pol2_fakerepeat/%s' ,conditions{z}));
    a = dir;

    bundle = 1000;
    [m,~] = size(a);
    j = 1;
    names = cell(1,m);
    for i = 1:m
        if ~isempty(strfind(a(i).name,thisfile))
            fprintf('On %s\n',conditions{z})
            fid = fopen(a(i).name,'r');
            if j == 1
                [s,w] = system(sprintf('wc -l %s',a(i).name));
                [len,~] = strread(w,'%d%s');
                A = zeros(m,len);
                C = cell(m,len);
            end
            k = 1;
            while ~feof(fid)
                B = textscan(fid,'%s%f',bundle,'HeaderLines',1);
                if k + bundle - 1 > len
                    A(j,k:k + length(B{2})-1) = B{2};
                    C(j,k:k+length(B{2})-1) = B{1};
                else
                    A(j,k:k+bundle-1) = B{2};
                    C(j,k:k+bundle-1) = B{1};
                    k = k + bundle;
                end   
            end
            names{j} = a(i).name;
            j = j + 1;
            fclose(fid);clear fid B
        end
    end
    Probs = A(1:j-1,:);
    readnames = C(1:j-1,:);
    names = names(1:j-1);
end