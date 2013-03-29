%Tell me if the peaks are within the gene body

A = eval(sprintf('%s_relativeto_%s',strcat(headers{headersind(2:end)}),headers{headersind(1)}));
subset = [A{1,2} A{1,3}];
chrsubset = A{1,1};

[m,n] = size(subset);

InGene = zeros(length(knownGene),m);
for i = 1:m
    left = subset(i,1);
    right = subset(i,2);
    knownGenesubset = find(strcmp(chrsubset(i),knownGene(:,4)));
    IND = 1;
    for j = 1:length(knownGenesubset)
        LeftEdge = cell2mat(knownGene(knownGenesubset(j),6))-1000;
        RightEdge = cell2mat(knownGene(knownGenesubset(j),7))+1000;
        win = (left >= LeftEdge & right <= RightEdge) | ...       %IP2 window fits inside of OS window
                (left <= LeftEdge & right >= LeftEdge) | ...         %IP2 window straddles start of OS window
                (left <= RightEdge & right >= RightEdge) | ...         %IP2 window staddles the end of the window
                (left <= LeftEdge & right >= RightEdge);              %OS window fits inside of the IP2 window
        if win == 1
           InGene(IND,i) = 1;
        else
           InGene(IND,i) = 0;
        end
    end
end