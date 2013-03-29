%Peak Percentage for different IPs that covers each gene
person = 'TY';
peak = 'Macs';

Table = cell(length(headers)+1,length(knownGene));

for j = 1:length(headers)
    x = importdata(sprintf('%s_%s_%s.txt',person,headers{j},peak));
    
    Text = x.textdata;
    Data = x.data;
    clear x
    
    for i = 1:length(knownGene)    
        Table(1,i) = knownGene{i,2};
        
        if cell2mat(knownGene(i,5))
            TSS = cell2mat(knownGene(i,6));
            TES = cell2mat(knownGene(i,7));
        else
            TSS = cell2mat(knownGene(i,7));
            TES = cell2mat(knownGene(i,6));
        end
        
        if strcmp(peak,'Lin_Conden')
        elseif strcmp(peak,'Poisson_Conden')
        else
        end
    end
end