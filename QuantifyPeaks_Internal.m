%Quantify the number of peaks that are on a gene/outside of a gene

%5' is -10kb-(-2kb), promoter is -2kb-1kb, gene body is 1kb-TES, 3' is
%TES-TES+10kb
control = find(strcmp('PI',headers));
chr = fieldnames(eval(headers{1}));

for j = 1:control-1
    
    %Initialize count values for each IP
    fiveprime = 0;
    promoter = 0;
    genebody = 0;
    threeprime = 0;
    outofgene = 0;
    total = 0;
    
    x = eval(sprintf('%s_Poiss_Conden',headers{j}));
    yy = eval(sprintf('%s_Chr_Conden',headers{j}));
    z = eval(sprintf('%s_TSS_Conden',headers{j}));
    
    for k = 1:length(chr)
        c = yy == k;
        cc = find(c);
        y = unique(x(cc(:),1:2),'rows');
        l = y(:,1);
        r = y(:,2);
        v  = zeros(1,length(l));
        total = total + length(l);
        
        f = unique(z(:,cc)','rows')';
        for i = 1:length(f)
            
            if f(1,i) - f(2,i) < 0
                
                %Calculate 5'
                a = l >= (f(1,i)-10000) & r < (f(1,i)-2000);
                a1 = r >= (f(1,i)-10000) & l < (f(1,i)-10000);
                a2 = l < (f(1,i)-2000) & r > (f(1,i)-2000);
                pr = abs(r(a1) - (f(1,i)-10000))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (f(1,i)-2000))./abs(r(a2)-l(a2));
                fiveprime = fiveprime + sum(a) + sum(pr) + sum(pl);
                v(a) = v(a) + 1;
                v(a1) = v(a1) + 1;
                v(a2) = v(a2) + 1;
                
                %Calculate promoter
                a = l >= (f(1,i)-2000) & r < (f(1,i)+1000);
                a1 = r >= (f(1,i)-2000) & l < (f(1,i)-2000);
                a2 = l < (f(1,i)+1000) & r > (f(1,i)+1000);
                pr = abs(r(a1) - (f(1,i)-2000))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (f(1,i)+1000))./abs(r(a2)-l(a2));
                promoter = promoter + sum(a) + sum(pr) + sum(pl);
                v(a) = v(a) + 1;
                v(a1) = v(a1) + 1;
                v(a2) = v(a2) + 1;
                
                %Calculate gene body
                a = l >= (f(1,i)+1000) & r < f(2,i);
                a1 = r >= (f(1,i)+1000) & l < (f(1,i)+1000);
                a2 = l < (f(2,i)) & r > f(2,i);
                pr = abs(r(a1) - (f(1,i)+1000))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (f(2,i)))./abs(r(a2)-l(a2));
                genebody = genebody +sum(a) + sum(pr) + sum(pl);
                v(a) = v(a) + 1;
                v(a1) = v(a1) + 1;
                v(a2) = v(a2) + 1;
                
                %Calculate 3'
                a = l >= (f(2,i)) & r < (f(2,i)+10000);
                a1 = r >= (f(2,i)) & l < f(2,i);
                a2 = l < (f(2,i)+10000) & r > (f(2,i)+10000);
                pr = abs(r(a1) - (f(2,i)))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (f(2,i)+10000))./abs(r(a2)-l(a2));
                threeprime = threeprime + sum(a) + sum(pr) + sum(pl);
                v(a) = v(a) + 1;
                v(a1) = v(a1) + 1;
                v(a2) = v(a2) + 1;
                
            else
                %Calculate 5'
                a = r <= (f(2,i)+10000) & l > (f(2,i)+2000);
                a1 = l <= (f(2,i)+10000) & r > (f(2,i)+10000);
                a2 = l < (f(2,i)+2000) & r  > (f(2,i)+2000);
                pl = abs(l(a1) - (f(2,i)+10000))./abs(r(a1)-l(a1));
                pr = abs(r(a2) - (f(2,i)+2000))./abs(r(a2)-l(a2));
                fiveprime = fiveprime + sum(a) + sum(pr) + sum(pl);
                v(a) = v(a) + 1;
                v(a1) = v(a1) + 1;
                v(a2) = v(a2) + 1;
                
                %Calculate promoter
                a = r <= (f(2,i)+2000) & l > (f(2,i)-1000);
                a1 = l <= (f(2,i)+2000) & r > (f(2,i)+2000);
                a2 = l < (f(2,i)-1000) & r  > (f(2,i)-1000);
                pl = abs(l(a1) - (f(2,i)+2000))./abs(r(a1)-l(a1));
                pr = abs(r(a2) - (f(2,i)-1000))./abs(r(a2)-l(a2));
                promoter = promoter + sum(a) + sum(pr) + sum(pl);
                v(a) = v(a) + 1;
                v(a1) = v(a1) + 1;
                v(a2) = v(a2) + 1;
                
                %Calculate gene body
                a = r <= (f(2,i)-1000) & l > f(1,i);
                a1 = l <= (f(2,i)-1000) & r > (f(2,i)-1000);
                a2 = l < (f(1,i)) & r  > f(1,i);
                pl = abs(l(a1) - (f(2,i)-1000))./abs(r(a1)-l(a1));
                pr = abs(r(a2) - (f(1,i)))./abs(r(a2)-l(a2));
                genebody = genebody +sum(a) + sum(pr) + sum(pl);
                v(a) = v(a) + 1;
                v(a1) = v(a1) + 1;
                v(a2) = v(a2) + 1;
                
                %Calculate 3'
                a = r <= (f(1,i)) & l > (f(1,i)-10000);
                a1 = l <= (f(1,i)) & r > f(1,i);
                a2 = l < (f(1,i)-10000) & r  > (f(1,i)-10000);
                pl = abs(l(a1) - (f(1,i)))./abs(r(a1)-l(a1));
                pr = abs(r(a2) - (f(1,i)-10000))./abs(r(a2)-l(a2));
                threeprime = threeprime + sum(a) + sum(pr) + sum(pl);
                v(a) = v(a) + 1;
                v(a1) = v(a1) + 1;
                v(a2) = v(a2) + 1;
            end
        end
        outofgene = outofgene + sum(v == 0);
    end
    figure(j)
    pie([fiveprime promoter genebody threeprime outofgene])
    title(sprintf('%s - Internal Norm.',headers{j}))
    legend('5 Prime','Promoter','Gene Body','3 Prime','Out of Gene')
    set(gca,'FontSize',13,'FontWeight','bold')
    text(-2,1,sprintf('Total Peaks = %s\n5 Prime = %s\nPromoter = %s\nGene Body = %s\n3 Prime = %s\nOut of Gene = %s',num2str(total),num2str(fiveprime),...
        num2str(promoter),num2str(genebody),num2str(threeprime),num2str(outofgene)),'FontSize',13','FontWeight','bold')
end

clear l r j k i pl pr a1 a2 a v x outofgene threeprime fiveprime genebody ...
    promoter total
