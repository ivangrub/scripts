function annotate_overlap(chr,A,per,ip,percent,overlap)

load ~/Documents/MATLAB/Tjian' Lab'/ChIP-Seq/ChIP-Seq' Scripts'/knownGene.mm9.mat
known = fopen('NewKnownGene.txt','r');
GI = textscan(known,'%s%s%d%s','delimiter','\t');
fclose(known);clear known

fid2 = fopen(sprintf('~/ChIP-Data/%s/Overlapping_%s_%s_%dPercentile.txt',per,overlap,ip,percent),'w');

for j = 1:length(A(:,1))
    Mid = mean([A(j,1) A(j,2)]);
   
    Avail = strcmp(knownGene(:,4),strcat(chr(j,1:5)));
    indexed = find(Avail);
    
    TSS = cell2mat(knownGene(Avail,5)) == 1;
    TSS_index = find(TSS);
    TES = cell2mat(knownGene(Avail,5)) == -1;
    TES_index = find(TES);
    
    up_tss_plus = Mid > cell2mat(knownGene(indexed(TSS_index),6));
    down_tss_plus = Mid < cell2mat(knownGene(indexed(TSS_index),6));
    up_tes_plus =  Mid > cell2mat(knownGene(indexed(TSS_index),7));
    down_tes_plus = Mid < cell2mat(knownGene(indexed(TSS_index),7));
    
    up_tss_minus = Mid > cell2mat(knownGene(indexed(TES_index),7));
    down_tss_minus = Mid < cell2mat(knownGene(indexed(TES_index),7));
    up_tes_minus =  Mid > cell2mat(knownGene(indexed(TES_index),6));
    down_tes_minus = Mid < cell2mat(knownGene(indexed(TES_index),6));
    
    
    if sum(up_tss_plus) > 0 || sum(up_tss_minus) > 0
        [dist_plus,k_plus] = min(abs(Mid - cell2mat(knownGene(indexed(TSS_index(up_tss_plus)),6))));
        [dist_minus,k_minus] = min(abs(Mid - cell2mat(knownGene(indexed(TES_index(up_tss_minus)),7))));
        if ~isempty(k_plus) && ~isempty(k_minus)
            [dist,ind] = min([dist_plus dist_minus]);
        elseif isempty(k_minus)
            dist = dist_plus;
            ind = 1;
        else
            dist = dist_minus;
            ind = 2;
        end
        if ind == 1
            f = find(up_tss_plus);
            k = indexed(TSS_index(f(k_plus)));
            strand = '+';
        else
            f = find(up_tss_minus);
            k = indexed(TES_index(f(k_minus)));
            strand = '-';
        end
        
        if ~isempty(k)
            a = strcmp(knownGene{k,1},GI{1,1});
            fprintf(fid2,'%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',char(knownGene{k,1}),char(knownGene{k,2}),GI{1,3}(a),char(knownGene{k,4}),strand,cell2mat(knownGene(k,6)),cell2mat(knownGene(k,7)),...
                Mid,dist,A(j,1),A(j,2),1,0,1,0);
        end
    end
    
    if sum(down_tss_plus) > 0 || sum(down_tss_minus) > 0
        [dist_plus,k_plus] = min(abs(Mid - cell2mat(knownGene(indexed(TSS_index(down_tss_plus)),6))));
        [dist_minus,k_minus] = min(abs(Mid - cell2mat(knownGene(indexed(TES_index(down_tss_minus)),7))));
        if ~isempty(k_plus) && ~isempty(k_minus)
            [dist,ind] = min([dist_plus dist_minus]);
        elseif isempty(k_minus)
            dist = dist_plus;
            ind = 1;
        else
            dist = dist_minus;
            ind = 2;
        end
        if ind == 1
            f = find(down_tss_plus);
            k = indexed(TSS_index(f(k_plus)));
            strand = '+';
        else
            f = find(down_tss_minus);
            k = indexed(TES_index(f(k_minus)));
            strand = '-';
        end
        if ~isempty(k)
            a = strcmp(knownGene{k,1},GI{1,1});
            fprintf(fid2,'%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',char(knownGene{k,1}),char(knownGene{k,2}),GI{1,3}(a),char(knownGene{k,4}),strand,cell2mat(knownGene(k,6)),cell2mat(knownGene(k,7)),...
                Mid,dist,A(j,1),A(j,2),1,0,0,1);
        end
    end
    if sum(up_tes_plus) > 0 || sum(up_tes_minus) > 0
        [dist_plus,k_plus] = min(abs(Mid - cell2mat(knownGene(indexed(TSS_index(up_tes_plus)),7))));
        [dist_minus,k_minus] = min(abs(Mid - cell2mat(knownGene(indexed(TES_index(up_tes_minus)),6))));
        if ~isempty(k_plus) && ~isempty(k_minus)
            [dist,ind] = min([dist_plus dist_minus]);
        elseif isempty(k_minus)
            dist = dist_plus;
            ind = 1;
        else
            dist = dist_minus;
            ind = 2;
        end
        if ind == 1
            f = find(up_tes_plus);
            k = indexed(TSS_index(f(k_plus)));
            strand = '+';
        else
            f = find(up_tes_minus);
            k = indexed(TES_index(f(k_minus)));
            strand = '-';
        end
        if ~isempty(k)
            a = strcmp(knownGene{k,1},GI{1,1});
            fprintf(fid2,'%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',char(knownGene{k,1}),char(knownGene{k,2}),GI{1,3}(a),char(knownGene{k,4}),strand,cell2mat(knownGene(k,6)),cell2mat(knownGene(k,7)),...
                Mid,dist,A(j,1),A(j,2),0,1,1,0);
        end
    end
    if sum(down_tes_plus) > 0 || sum(down_tes_minus) > 0
        [dist_plus,k_plus] = min(abs(Mid - cell2mat(knownGene(indexed(TSS_index(down_tes_plus)),7))));
        [dist_minus,k_minus] = min(abs(Mid - cell2mat(knownGene(indexed(TES_index(down_tes_minus)),6))));
        if ~isempty(k_plus) && ~isempty(k_minus)
            [dist,ind] = min([dist_plus dist_minus]);
        elseif isempty(k_minus)
            dist = dist_plus;
            ind = 1;
        else
            dist = dist_minus;
            ind = 2;
        end
        if ind == 1
            f = find(down_tes_plus);
            k = indexed(TSS_index(f(k_plus)));
            strand = '+';
        else
            f = find(down_tes_minus);
            k = indexed(TES_index(f(k_minus)));
            strand = '-';
        end
        
        if ~isempty(k)
            a = strcmp(knownGene{k,1},GI{1,1});
            fprintf(fid2,'%s\t%s\t%d\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',char(knownGene{k,1}),char(knownGene{k,2}),GI{1,3}(a),char(knownGene{k,4}),strand,cell2mat(knownGene(k,6)),cell2mat(knownGene(k,7)),...
                Mid,dist,A(j,1),A(j,2),0,1,0,1);
    
        end
    end
end
fclose(fid2);clear fid2