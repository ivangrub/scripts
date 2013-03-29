%Quantify the number of peaks that are on a gene/outside of a gene

%5' is -10kb-(-2kb), promoter is -2kb-1kb, gene body is 1kb-TES, 3' is TES-TES+10kb

control = find(strcmp('PI',headers));
chr = fieldnames(eval(headers{1}));

for j = 1:control-1
	
	%Initialize count values for each IP
	fiveprime = 0;
	promoter = 0;
	genebody = 0;
	threeprime = 0;
	outofgene = 0;
    
	x = eval(sprintf('%s_Poiss_Conden',headers{j}));
    
	for k = 1:length(chr)
		c = Chr_conden == k;
		cc = find(c);
        y = unique(x(cc(:),1:2),'rows');
        l = y(:,1);
        r = y(:,2);
        
        f = unique(TSS_conden(:,cc)','rows')';
		for i = 1:length(f)
            
            if f(1,i) - f(2,i) < 0
                
                %Calculate 5' 
                a = l >= (f(1,i)-10000) & r < (f(1,i)-2000);
                a1 = r >= (f(1,i)-10000) & l < (f(1,i)-10000);
                a2 = l < (f(1,i)-2000) & r > (f(1,i)-2000);
                pr = abs(r(a1) - (f(1,i)-10000))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (f(1,i)-2000))./abs(r(a2)-l(a2));
                fiveprime = fiveprime + sum(a) + sum(pr) + sum(pl);

                %Calculate promoter
                a = l >= (f(1,i)-2000) & r < (f(1,i)+1000);
                a1 = r >= (f(1,i)-2000) & l < (f(1,i)-2000);
                a2 = l < (f(1,i)+1000) & r > (f(1,i)+1000);
                pr = abs(r(a1) - (f(1,i)-2000))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (f(1,i)+1000))./abs(r(a2)-l(a2));
                promoter = promoter + sum(a) + sum(pr) + sum(pl);

                %Calculate gene body
                a = l >= (f(1,i)+1000) & r < f(2,i);
                a1 = r >= (f(1,i)+1000) & l < (f(1,i)+1000);
                a2 = l < (f(2,i)) & r > f(2,i);
                pr = abs(r(a1) - (f(1,i)+1000))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (f(2,i)))./abs(r(a2)-l(a2));
                genebody = genebody +sum(a) + sum(pr) + sum(pl);
                
                %Calculate 3'
                a = l >= (f(2,i)) & r < (f(2,i)+10000);
                a1 = r >= (f(2,i)) & l < f(2,i);
                a2 = l < (f(2,i)+10000) & r > (f(2,i)+10000);
                pr = abs(r(a1) - (f(2,i)))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (f(2,i)+10000))./abs(r(a2)-l(a2));
                threeprime = threeprime + sum(a) + sum(pr) + sum(pl);
                
                %Calculate out of gene
                a = r < f(1,i)-10000 | l > (f(2,i)+10000);
                outofgene = outofgene + sum(a);
            else
                %Calculate 5' 
                a = r <= (f(2,i)+10000) & l > (f(2,i)+2000);
                a1 = l <= (f(2,i)+10000) & r > (f(2,i)+10000);
                a2 = l < (f(2,i)+2000) & r  > (f(2,i)+2000);
                pr = abs(l(a1) - (f(2,i)+10000))./abs(r(a1)-l(a1));
                pl = abs(r(a2) - (f(2,i)+2000))./abs(r(a2)-l(a2));
                fiveprime = fiveprime + sum(a) + sum(pr) + sum(pl);

                %Calculate promoter
                a = r <= (f(2,i)+2000) & l > (f(2,i)-1000);
                a1 = l <= (f(2,i)+2000) & r > (f(2,i)+2000);
                a2 = l < (f(2,i)-1000) & r  > (f(2,i)-1000);
                pr = abs(l(a1) - (f(2,i)+2000))./abs(r(a1)-l(a1));
                pl = abs(r(a2) - (f(2,i)-1000))./abs(r(a2)-l(a2));
                promoter = promoter + sum(a) + sum(pr) + sum(pl);

                %Calculate gene body
                a = r <= (f(2,i)-1000) & l > f(1,i);
                a1 = l <= (f(2,i)-1000) & r > (f(2,i)-1000);
                a2 = l < (f(1,i)) & r  > f(1,i);
                pr = abs(l(a1) - (f(2,i)-1000))./abs(r(a1)-l(a1));
                pl = abs(r(a2) - (f(1,i)))./abs(r(a2)-l(a2));
                genebody = genebody +sum(a) + sum(pr) + sum(pl);
                
                %Calculate 3'
                a = r <= (f(1,i)) & l > (f(1,i)-10000);
                a1 = l <= (f(1,i)) & r > f(1,i);
                a2 = l < (f(1,i)-10000) & r  > (f(1,i)-10000);
                pr = abs(l(a1) - (f(1,i)))./abs(r(a1)-l(a1));
                pl = abs(r(a2) - (f(1,i)-10000))./abs(r(a2)-l(a2));
                threeprime = threeprime + sum(a) + sum(pr) + sum(pl);
                
                %Calculate out of gene
                a = r < f(1,i)-10000 | l > (f(2,i)+10000);
                outofgene = outofgene + sum(a);
            end
		end
	end
end