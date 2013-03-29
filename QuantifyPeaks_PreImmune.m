%Quantify the number of peaks that are on a gene/outside of a gene

MaxGenes = 1.1;           %Does not include this value
NumofSections = 1;
WOutofgene = 0;


%5' is -10kb-(-2kb), promoter is -2kb-1kb, gene body is 1kb-TES, 3' is
%TES-TES+10kb

chr = fieldnames(eval(sprintf('%s.win',headers{IPidx(1)})));
TSS = cell(length(knownGene),4);
for i = 1:length(knownGene)
    if cell2mat(knownGene(i,5)) == 1
        TSS(i,:) = {knownGene{i,4} cell2mat(knownGene(i,5)) cell2mat(knownGene(i,6)) cell2mat(knownGene(i,7))};
    else
        TSS(i,:) = {knownGene{i,4} cell2mat(knownGene(i,5)) cell2mat(knownGene(i,7)) cell2mat(knownGene(i,6))};
    end
end

for j = 1:length(IPidx)
    
    %Initialize count values for each IP
    Count = zeros(1,4);
    
    x = eval(sprintf('%s_condensed',headers{IPidx(j)}));
    yy = eval(sprintf('%s_Chr',headers{IPidx(j)}));
    z = eval(sprintf('%s_TSS',headers{IPidx(j)}));
    fullmat = [yy' x];
    IndivPeaks = unique(fullmat(:,1:3),'rows');
    total = length(IndivPeaks);
    
    for k = 1:length(chr)
        c = IndivPeaks(:,1) == k;
        cc = find(c);
        l = IndivPeaks(cc,2);
        r = IndivPeaks(cc,3);
        vf_fp  = zeros(1,length(l));
        vf_pr  = zeros(1,length(l));
        vf_gb  = zeros(1,length(l));
        vf_tp  = zeros(1,length(l));
        vr_fp = zeros(1,length(l));
        vr_pr = zeros(1,length(l));
        vr_gb = zeros(1,length(l));
        vr_tp = zeros(1,length(l));
        
        
        ind = find(strcmp(chr{k},TSS(:,1)));
        f = cell2mat(TSS(ind,2:4))';
        if ~isempty(f)
            for i = 1:length(f(1,:))
                
                if f(1,i) == 1
                    
                    %Calculate 5'
                    a = l >= (f(2,i)-10000) & r < (f(2,i)-2000);
                    a1 = r >= (f(2,i)-10000) & l < (f(2,i)-10000);
                    a2 = l < (f(2,i)-2000) & r > (f(2,i)-2000);
                    pr = abs(r(a1) - (f(2,i)-10000))./abs(r(a1)-l(a1));
                    pl = abs(l(a2) - (f(2,i)-2000))./abs(r(a2)-l(a2));
                    
                    vf_fp(a) = vf_fp(a) + 1;
                    if ~isempty(pr)
                        a = find(a1);
                        for z = 1:length(a)
                            vf_fp(a(z)) = vf_fp(a(z)) + pr(z);
                        end
                    end
                    if ~isempty(pl)
                        a = find(a2);
                        for z = 1:length(a)
                            vf_fp(a(z)) = vf_fp(a(z)) + pl(z);
                        end
                    end
                    
                    %Calculate promoter
                    a = l >= (f(2,i)-2000) & r < (f(2,i)+1000);
                    a1 = r >= (f(2,i)-2000) & l < (f(2,i)-2000);
                    a2 = l < (f(2,i)+1000) & r > (f(2,i)+1000);
                    pr = abs(r(a1) - (f(2,i)-2000))./abs(r(a1)-l(a1));
                    pl = abs(l(a2) - (f(2,i)+1000))./abs(r(a2)-l(a2));
                    
                    vf_pr(a) = vf_pr(a) + 1;
                    if ~isempty(pr)
                        a = find(a1);
                        for z = 1:length(a)
                            vf_pr(a(z)) = vf_pr(a(z)) + pr(z);
                        end
                    end
                    if ~isempty(pl)
                        a = find(a2);
                        for z = 1:length(a)
                            vf_pr(a(z)) = vf_pr(a(z)) + pl(z);
                        end
                    end
                    
                    %Calculate gene body
                    a = l >= (f(2,i)+1000) & r < f(3,i);
                    a1 = r >= (f(2,i)+1000) & l < (f(2,i)+1000);
                    a2 = l < (f(3,i)) & r > f(3,i);
                    pr = abs(r(a1) - (f(2,i)+1000))./abs(r(a1)-l(a1));
                    pl = abs(l(a2) - (f(3,i)))./abs(r(a2)-l(a2));
                    
                    vf_gb(a) = vf_gb(a) + 1;
                    if ~isempty(pr)
                        a = find(a1);
                        for z = 1:length(a)
                            vf_gb(a(z)) = vf_gb(a(z)) + pr(z);
                        end
                    end
                    if ~isempty(pl)
                        a = find(a1);
                        for z = 1:length(a)
                            vf_gb(a(z)) = vf_gb(a(z)) + pl(z);
                        end
                    end
                    
                    %Calculate 3'
                    a = l >= (f(3,i)) & r < (f(3,i)+10000);
                    a1 = r >= (f(3,i)) & l < f(3,i);
                    a2 = l < (f(3,i)+10000) & r > (f(3,i)+10000);
                    pr = abs(r(a1) - (f(3,i)))./abs(r(a1)-l(a1));
                    pl = abs(l(a2) - (f(3,i)+10000))./abs(r(a2)-l(a2));
                    
                    vf_tp(a) = vf_tp(a) + 1;
                    if ~isempty(pr)
                        a = find(a1);
                        for z = 1:length(a)
                            vf_tp(a(z)) = vf_tp(a(z)) + pr(z);
                        end
                    end
                    if ~isempty(pl)
                        a = find(a2);
                        for z = 1:length(a)
                            vf_tp(a(z)) = vf_tp(a(z)) + pl(z);
                        end
                    end
                else
                    %Calculate 5'
                    a = r <= (f(2,i)+10000) & l > (f(2,i)+2000);
                    a1 = l <= (f(2,i)+10000) & r > (f(2,i)+10000);
                    a2 = l < (f(2,i)+2000) & r  > (f(2,i)+2000);
                    pl = abs(l(a1) - (f(2,i)+10000))./abs(r(a1)-l(a1));
                    pr = abs(r(a2) - (f(2,i)+2000))./abs(r(a2)-l(a2));
                    
                    vr_fp(a) = vr_fp(a) + 1;
                    if ~isempty(pl)
                        a = find(a1);
                        for z = 1:length(a)
                            vr_fp(a(z)) = vr_fp(a(z)) + pl(z);
                        end
                    end
                    if ~isempty(pr)
                        a = find(a2);
                        for z = 1:length(a)
                            vr_fp(a(z)) = vr_fp(a(z)) + pr(z);
                        end
                    end
                    
                    %Calculate promoter
                    a = r <= (f(2,i)+2000) & l > (f(2,i)-1000);
                    a1 = l <= (f(2,i)+2000) & r > (f(2,i)+2000);
                    a2 = l < (f(2,i)-1000) & r  > (f(2,i)-1000);
                    pl = abs(l(a1) - (f(2,i)+2000))./abs(r(a1)-l(a1));
                    pr = abs(r(a2) - (f(2,i)-1000))./abs(r(a2)-l(a2));
                    
                    vr_pr(a) = vr_pr(a) + 1;
                    if ~isempty(pl)
                        a = find(a1);
                        for z = 1:length(a)
                            vr_pr(a(z))= vr_pr(a(z)) + pl(z);
                        end
                    end
                    if ~isempty(pr)
                        a = find(a2);
                        for z = 1:length(a)
                            vr_pr(a(z)) = vr_pr(a(z)) + pr(z);
                        end
                    end
                    
                    %Calculate gene body
                    a = r <= (f(2,i)-1000) & l > f(3,i);
                    a1 = l <= (f(2,i)-1000) & r > (f(2,i)-1000);
                    a2 = l < (f(3,i)) & r  > f(3,i);
                    pl = abs(l(a1) - (f(2,i)-1000))./abs(r(a1)-l(a1));
                    pr = abs(r(a2) - (f(3,i)))./abs(r(a2)-l(a2));
                    
                    vr_gb(a) = vr_gb(a) + 1;
                    if ~isempty(pl)
                        a = find(a1);
                        for z = 1:length(a)
                            vr_gb(a(z)) = vr_gb(a(z)) + pl(z);
                        end
                    end
                    if ~isempty(pr)
                        a = find(a2);
                        for z = 1:length(a)
                            vr_gb(a(z)) = vr_gb(a(z)) + pr(z);
                        end
                    end
                    
                    %Calculate 3'
                    a = r <= (f(3,i)) & l > (f(3,i)-10000);
                    a1 = l <= (f(3,i)) & r > f(3,i);
                    a2 = l < (f(3,i)-10000) & r  > (f(3,i)-10000);
                    pl = abs(l(a1) - (f(3,i)))./abs(r(a1)-l(a1));
                    pr = abs(r(a2) - (f(3,i)-10000))./abs(r(a2)-l(a2));
                    
                    vr_tp(a) = vr_tp(a) + 1;
                    if ~isempty(pl)
                        a = find(a1);
                        for z = 1:length(a)
                            vr_tp(a(z)) = vr_tp(a(z)) + pl(z);
                        end
                    end
                    if ~isempty(pr)
                        a = find(a2);
                        for z = 1:length(a)
                            vr_tp(a(z)) = vr_tp(a(z)) + pr(z);
                        end
                    end
                end
            end
        end
        vf = [vf_fp; vf_pr ;vf_gb; vf_tp];
        vr = [vr_fp; vr_pr ;vr_gb; vr_tp];
        af = vf == 0;
        bf = sum(af);
        xf = bf==4;
        ar = vr(:,xf) == 0;
        br = sum(ar);
        xr = br == 4;
        outofgene = sum(xr);
        for i = 1:length(vr(1,:));
            a = (vf(:,i) > 0 & vf(:,i) <= MaxGenes) | (vr(:,i) > 0 & vr(:,i) <= MaxGenes);
            if sum(a) <= NumofSections && sum(a) > 0
                idx = find(a);
                Count(1,idx) = Count(1,idx) + sum(vr(idx,i)) + sum(vf(idx,i));
            end
        end
    end
    figure(j)
    
    if WOutofgene == 1
        pie([Count(1) Count(2) Count(3) Count(4) outofgene])
        title(sprintf('%s - Pre-Immune',strrep(headers{IPidx(j)},'_',' ')))
        legend('5 Prime','Promoter','Gene Body','3 Prime','Out of Gene')
        set(gca,'FontSize',13,'FontWeight','bold')
        text(-2,1,sprintf('Total Peaks = %s\n5 Prime = %s\nPromoter = %s\nGene Body = %s\n3 Prime = %s\nOut of Gene = %s',num2str(total),num2str(Count(1)),...
            num2str(Count(2)),num2str(Count(3)),num2str(Count(4)),num2str(outofgene)),'FontSize',13','FontWeight','bold')
    else
        pie([Count(1) Count(2) Count(3) Count(4)])
        title(sprintf('%s - Pre-Immune',strrep(headers{IPidx(j)},'_',' ')))
        legend('5 Prime','Promoter','Gene Body','3 Prime')
        set(gca,'FontSize',13,'FontWeight','bold')
        text(-2,1,sprintf('Total Peaks = %s\n5 Prime = %s\nPromoter = %s\nGene Body = %s\n3 Prime = %s',num2str(total),num2str(Count(1)),...
            num2str(Count(2)),num2str(Count(3)),num2str(Count(4))),'FontSize',13','FontWeight','bold')
    end
end

clear l r j k i pl pr a1 a2 v x outofgene threeprime fiveprime genebody ...
    promoter total

