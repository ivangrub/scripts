%Compare ES and MEF cell peak overlaps

clc

merge = 1;

IP1 = PolII_N_relativeto_TAF7;
IP2 = MEF_PolII_relativeto_MEF_TAF7;

IND = 1;
for k = 1:length(IP1{1,1})
    chr = strcmp(IP2{1,1}',IP1{1,1}(k));
    win = (IP2{1,2} >= IP1{1,2}(k)-width & IP2{1,3} <= IP1{1,3}(k) + width) | ...       %IP2 window fits inside of IP1 window
        (IP2{1,2} < IP1{1,2}(k)-width & IP2{1,3} >= IP1{1,2}(k) - width) | ...         %IP2 window straddles start of IP1 window
        (IP2{1,2} <= IP1{1,3}(k)+width & IP2{1,3} > IP1{1,3}(k) + width) | ...         %IP2 window staddles the end of IP1 window
        (IP2{1,2} <= IP1{1,2}(k)-width & IP2{1,3} >= IP1{1,3}(k) + width);              %IP1 window fits inside of the IP2 window
    Overlap = (chr & win);
    if sum(Overlap) == 1
        chrsubset(IND) = IP1{1,1}(k);
        subset(IND,[1 2 3 4]) = [IP1{1,2}(k) IP1{1,3}(k) IP2{1,2}(chr & win) IP2{1,3}(chr & win)];
        IND = IND + 1;
        continue
    elseif sum(Overlap) > 1
        
        chrsubset(IND:IND+sum(Overlap)-1) = IP1{1,1}(k);
        
        subset(IND:IND+sum(Overlap)-1,1) = IP1{1,2}(k);
        subset(IND:IND+sum(Overlap)-1,2) = IP1{1,3}(k);
        
        
        subset(IND:IND+sum(Overlap)-1,[3 4]) = [IP2{1,2}(chr&win) IP2{1,3}(chr&win)];
        IND = IND + sum(Overlap);
        
        continue
    elseif sum(Overlap) == 0
        chrsubset(IND) = IP1{1,1}(k);
        subset(IND,[1 2 3 4]) = [IP1{1,2}(k) IP1{1,3}(k) 0 0];
        IND = IND + 1;
        continue
    end
end

cover = subset(:,3) ~= 0;
nocover = subset(:,3) == 0;

same = subset(cover,:);
same_chr = chrsubset(cover);

diff = subset(nocover,:);
diff_chr = chrsubset(nocover);

clear subset

name = {'same' 'diff'};
for j = 1:2
    x = eval(name{j});
    win = zeros(length(x),4);
    C = zeros(length(x),2);
    
    y = eval(sprintf('%s_chr',name{j}));
    
    if merge == 0
        win(:,1) = (x(:,3) >= x(:,1) & x(:,4) <= x(:,2)); 
        win(:,2) = (x(:,3) < x(:,1) & x(:,4) >= x(:,1));
        win(:,3) = (x(:,3) <= x(:,2) & x(:,4) > x(:,2));
        win(:,4) = (x(:,3) <= x(:,1) & x(:,4) >= x(:,2));
        win = logical(win);
        C(win(:,1),1) = x(win(:,1),3);
        C(win(:,1),2) = x(win(:,1),4);
        C(win(:,2),1) = x(win(:,2),1);
        C(win(:,2),2) = x(win(:,2),4);
        C(win(:,3),1) = x(win(:,3),3);
        C(win(:,3),2) = x(win(:,3),2);
        C(win(:,4),1) = x(win(:,4),1);
        C(win(:,4),2) = x(win(:,4),2);
    else
        C(:,1) = min(x,[],2);
        C(:,2) = max(x,[],2);
    end
    
    if j == 2
        C = x(:,1:2);
    end
    
    clear win
    [m,n] = size(x);
    
    fid = fopen(sprintf('TY_ReltoESC_merge_2IP_%s.txt',name{j}),'w');
    InGene = zeros(length(knownGene),m);
    for i = 1:m
        knownGenesubset = strcmp(y(i),knownGene(:,4));
        left = C(i,1);
        right = C(i,2);
        LeftEdge = cell2mat(knownGene(:,6))-1000;
        RightEdge = cell2mat(knownGene(:,7))+1000;
        win = (left >= LeftEdge & right <= RightEdge) | ...       %IP2 window fits inside of OS window
            (left < LeftEdge & right >= LeftEdge) | ...         %IP2 window straddles start of OS window
            (left <= RightEdge & right > RightEdge) | ...         %IP2 window staddles the end of the window
            (left <= LeftEdge & right >= RightEdge);              %OS window fits inside of the IP2 window
        InGene(:,i) = win & knownGenesubset;
         if sum(InGene(:,i)) == 1
            fprintf(fid,'%s\t%s\t%s\t%d\t%d\n',knownGene{logical(InGene(:,i)),4},knownGene{logical(InGene(:,i)),2},knownGene{logical(InGene(:,i)),3},C(i,1),C(i,2));
        elseif sum(InGene(:,i)) > 1
            gs = find(InGene(:,i));
            for g = 1:length(gs)
                fprintf(fid,'%s\t%s\t%s\t%d\t%d\n',knownGene{gs(g),4},knownGene{gs(g),2},knownGene{gs(g),3},C(i,1),C(i,2));
             end
        elseif sum(InGene(:,i)) == 0
            fprintf(fid,'%s\t%s\t%s\t%d\t%d\n',y{i},'NaN','NaN',C(i,1),C(i,2));
        end
    end
    fclose(fid); clear fid InGene
    

end