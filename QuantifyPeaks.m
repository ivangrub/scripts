%Quantify the number of peaks that are on a gene/outside of a gene

%5' is -10kb-(-2kb), promoter is -2kb-1kb, gene body is 1kb-TES, 3' is TES-TES+10kb

control = find(strcmp('PI',headers));
IP = eval(headers{1});
chr = fieldnames(IP.win);

for j = 1%:control-1
	
	%Initialize count values for each IP
	fiveprime = 0;
	promoter = 0;
	genebody = 0;
	threeprime = 0;
	total = 0;
    
	x = eval(sprintf('%s_condensed',headers{j}));
    
	for k = 1:length(chr)
		c = find(strcmp(knownGene(:,4),chr{k}));
        f = unique(cell2mat(knownGene(c,6:7)),'rows');
		cc = find(Chr_conden == k);
        y = unique(x(cc(:),1:2),'rows');
        l = y(:,1);
        r = y(:,2);
        total = total + length(r);
		for i = 1:length(f)
            
            if cell2mat(knownGene(c(i),5)) == 1
                
                %Calculate 5' 
                a = l >= (knownGene{c(i),6}-10000) & r < (knownGene{c(i),6}-2000);
                a1 = r > (knownGene{c(i),6}-10000) & l < (knownGene{c(i),6}-10000);
                a2 = l < (knownGene{c(i),6}-2000) & r > (knownGene{c(i),6}-2000);
                pr = abs(r(a1) - (knownGene{c(i),6}-10000))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (knownGene{c(i),6}-2000))./abs(r(a2)-l(a2));
                fiveprime = fiveprime + sum(a) + sum(pr) + sum(pl);

                %Calculate promoter
                a = l >= (knownGene{c(i),6}-2000) & r < (knownGene{c(i),6}+1000);
                a1 = r > (knownGene{c(i),6}-2000) & l < (knownGene{c(i),6}-2000);
                a2 = l < (knownGene{c(i),6}+1000) & r > (knownGene{c(i),6}+1000);
                pr = abs(r(a1) - (knownGene{c(i),6}-2000))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (knownGene{c(i),6}+1000))./abs(r(a2)-l(a2));
                promoter = promoter + sum(a) + sum(pr) + sum(pl);

                %Calculate gene body
                a = l >= (knownGene{c(i),6}+1000) & r < knownGene{c(i),7};
                a1 = r > (knownGene{c(i),6}+1000) & l < knownGene{c(i),6}+1000;
                a2 = l < (knownGene{c(i),7}) & r > knownGene{c(i),7};
                pr = abs(r(a1) - (knownGene{c(i),6}+1000))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (knownGene{c(i),7}))./abs(r(a2)-l(a2));
                genebody = genebody +sum(a) + sum(pr) + sum(pl);
                
                %Calculate 3'
                a = l >= (knownGene{c(i),7}) & r < (knownGene{c(i),7}+10000);
                a1 = r > (knownGene{c(i),7}) & l < knownGene{c(i),7};
                a2 = l < (knownGene{c(i),7}+10000) & r > (knownGene{c(i),7}+10000);
                pr = abs(r(a1) - (knownGene{c(i),7}))./abs(r(a1)-l(a1));
                pl = abs(l(a2) - (knownGene{c(i),7}+10000))./abs(r(a2)-l(a2));
                threeprime = threeprime + sum(a) + sum(pr) + sum(pl);
            else
                %Calculate 5' 
                a = r <= (knownGene{c(i),7}+10000) & l > (knownGene{c(i),7}+2000);
                a1 = l < (knownGene{c(i),7}+10000) & r > knownGene{c(i),7}+10000;
                a2 = l < (knownGene{c(i),7}+2000) & r  > knownGene{c(i),7}+2000;
                pl = abs(l(a1) - (knownGene{c(i),7}+10000))./abs(r(a1)-l(a1));
                pr = abs(r(a2) - (knownGene{c(i),7}+2000))./abs(r(a2)-l(a2));
                fiveprime = fiveprime + sum(a) + sum(pr) + sum(pl);

                %Calculate promoter
                a = r <= (knownGene{c(i),7}+2000) & l > (knownGene{c(i),7}-1000);
                a1 = l < (knownGene{c(i),7}+2000) & r > (knownGene{c(i),7}+2000);
                a2 = l < (knownGene{c(i),7}-1000) & r  > (knownGene{c(i),7}-1000);
                pl = abs(l(a1) - (knownGene{c(i),7}+2000))./abs(r(a1)-l(a1));
                pr = abs(r(a2) - (knownGene{c(i),7}-1000))./abs(r(a2)-l(a2));
                promoter = promoter + sum(a) + sum(pr) + sum(pl);

                %Calculate gene body
                a = r <= (knownGene{c(i),7}-1000) & l > knownGene{c(i),6};
                a1 = l < (knownGene{c(i),7}-1000) & r > (knownGene{c(i),7}-1000);
                a2 = l < (knownGene{c(i),6}) & r  > knownGene{c(i),6};
                pl = abs(l(a1) - (knownGene{c(i),7}-1000))./abs(r(a1)-l(a1));
                pr = abs(r(a2) - (knownGene{c(i),6}))./abs(r(a2)-l(a2));
                genebody = genebody +sum(a) + sum(pr) + sum(pl);
                
                %Calculate 3'
                a = r <= (knownGene{c(i),6}) & l > (knownGene{c(i),6}-10000);
                a1 = l < (knownGene{c(i),6}) & r > knownGene{c(i),6};
                a2 = l < (knownGene{c(i),6}-10000) & r  > (knownGene{c(i),6}-10000);
                pl = abs(l(a1) - (knownGene{c(i),6}))./abs(r(a1)-l(a1));
                pr = abs(r(a2) - (knownGene{c(i),6}-10000))./abs(r(a2)-l(a2));
                threeprime = threeprime + sum(a) + sum(pr) + sum(pl);
            end
		end
	end
end