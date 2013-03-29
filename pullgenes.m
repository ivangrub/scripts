%Pull out the genes in each of the clusters

function pullgenes(group,kcluster,uniquegene,wcpggene,restgene)

n = max(kcluster);
if group == 1
    for i = 1:n
        m = kcluster == i;
        assignin('base',sprintf('kgene%d',i),uniquegene(m,1:2));
    end
elseif group == 2
    for i = 1:n
       m = kcluster == i;
       assignin('base',sprintf('kgene%d',i),wcpggene(m,1:2));
    end
else
    for i = 1:n
       m = kcluster == i;
       assignin('base',sprintf('kgene%d',i),restgene(m,1:2));
    end
end
end
