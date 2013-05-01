
conditions = {'pol2_100bp_neigh1_wB','pol2_300bp_wB','pol2_300bp_neigh1_wB'};
thisfile = 'ReadProbability.txt';

names = cell(1,length(conditions));
for z = 1:length(conditions)
    cd(sprintf('/Volumes/Genomic1/IG_express/fly_pol2_k100_fakerepeat/%s' ,conditions{z}));
    a = dir;

    bundle = 1000;
    [m,~] = size(a);
    j = 1;
    
    for i = 1:m
        if ~isempty(strfind(a(i).name,thisfile))
            fprintf('On %s\n',conditions{z})
            fid = fopen(a(i).name,'r');
            if z == 1
                [s,w] = system(sprintf('wc -l %s',a(i).name));
                [len,~] = strread(w,'%d%s');
                Probs = zeros(length(conditions),len);
                readnames = zeros(length(conditions),len);
                A = zeros(m,len);
                C = zeros(m,len);
            end
            k = 1;
            while ~feof(fid)
                B = textscan(fid,'%s%f',bundle,'HeaderLines',1);
                if k + bundle - 1 > len
                    A(j,k:k + length(B{2})-1) = B{2};
                    C(j,k:k+length(B{1})-1) = str2double(strrep(B{1},'SRR',''));
                else
                    A(j,k:k+bundle-1) = B{2};
                    C(j,k:k+bundle-1) = str2double(strrep(B{1},'SRR',''));
                    k = k + bundle;
                end   
            end
            j = j + 1;
            fclose(fid);clear fid B
        end
    end
    Probs(z,:) = A(1:j-1,:);
    readnames(z,:) = C(1:j-1,:);
    names(1,z) = conditions(z);
end