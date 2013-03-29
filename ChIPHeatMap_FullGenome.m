% Make heat maps of list of 5kb regions around a TSS for inputted genes

headers = {'ES_TAF7' 'ES_TBP' 'ES_PolII' 'ES_TAF1' 'MEF_TAF7' 'MEF_TBP' 'MEF_PolII'};

genes = fopen('refGene.mm9.txt');
glist = textscan(genes,'%s%s%s%d%d%d%d%d','delimiter','\t');
fclose(genes); clear genes

A = zeros(length(glist{1}),10000/window);
[m,n] = size(A);
chr = glist{2};


for i = 1:length(headers)
    x = eval(headers{i});
    k = 1;
    for j = 1:length(glist{1})
        if strcmp(glist{3}(j),'+')
            st = double(round(glist{4}(j)/window));
        else
            st = double(round(glist{5}(j)/window));
        end 
        if ~isempty(strfind(chr{j},'random'))
            continue
        end
        if st-n/2 < 0
            A(k,:) = x.win.(chr{j})(2,st:st+n-1);
        elseif st+n/2-1 > length(x.win.(chr{j})(2,:))
            add = st + n/2-1 - length(x.win.(chr{j})(2,:));
            A(k,:) = x.win.(chr{j})(2,st-(n/2+add):st+n/2-1-add);
        else
            A(k,:) = x.win.(chr{j})(2,st-n/2:st+n/2-1);
        end
        k = k + 1;
    end
    
    figure(i)
    imagesc(log2(A))
    title(strrep(headers{i},'_',' '))
    set(gca,'XTickLabel',-3750:1250:5000)
    saveas(gca,sprintf('%s_AllGenes.pdf',headers{i}))
    
end


