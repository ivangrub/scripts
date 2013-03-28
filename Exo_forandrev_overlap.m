offset = 0;

new = fopen('NewPeaks.bed','w');
coord = cell2mat(TFIIH(:,3:4));

chr = unique(TFIIH(:,1));
forward = TFIIH(strcmp(TFIIH(:,2),'+'),:);
fcoord = coord(strcmp(TFIIH(:,2),'+'),:);
reverse = TFIIH(strcmp(TFIIH(:,2),'-'),:);
freverse = coord(strcmp(TFIIH(:,2),'-'),:);

for i = 1:length(chr)
    ind = strcmp(chr(i,1),forward(:,1));
    ind2 = strcmp(chr(i,1),reverse(:,1));
    fp = fcoord(ind,1:2);
    rp = freverse(ind2,1:2);
    [n,~] = size(fp);
    for j = 1:n
        a = (fp(j,2) + offset >= rp(:,1) & fp(j,2) + offset <= rp(:,2));
        if sum(a) > 0
            ind = find(fp(j,2) >= rp(:,1) & fp(j,2) <= rp(:,2));
            for k = 1:length(ind)
                fprintf(new,'%s\t%d\t%d\n',char(chr{i,1}),fp(j,1),rp(ind(k),2));
            end
        end
    end
end

fclose(new);clear new