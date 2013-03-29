gene = 'Pou5f1';


g = find(strcmp(gene,knownGene(:,2)),1,'first');
chr = strcat('chr',knownGene{g,4});
x1 = cell2mat(knownGene(g,6))-5000;
x2 = cell2mat(knownGene(g,7))+5000;


win = floor(x1/25):floor(x2/25);

lin_c = find(strcmp(chr,Lin_Conden_text(1:end-1,5)));
linc = find(Lin_Conden(lin_c,4) > x1 & Lin_Conden(lin_c,5) < x2);
lin_s = find(strcmp(chr,Lin_Sig_text(1:end-1,5)));
lins = find(Lin_Sig(lin_s,3) > x1 & Lin_Sig(lin_s,4) < x2);
poi_c = find(strcmp(chr,Poisson_Conden_text(1:end-1,5)));
poic = find(Poisson_Conden(poi_c,4) > x1 & Poisson_Conden(poi_c,5) < x2);
poi_s = find(strcmp(chr,Poisson_Sig_text(1:end-1,5)));
pois = find(Poisson_Sig(poi_s,3) > x1 & Poisson_Sig(poi_s,4) < x2);
mac = find(strcmp(chr,Macs_text(1:end-1,1)));
macs = find(str2double(Macs_text(mac,2)) > x1-3000 & str2double(Macs_text(mac,3)) < x2);

figure(1)
subplot(2,1,1),hold on
if isempty(lins) == 0
    for i = 1:length(lins)
         sig = Lin_Sig(lin_s(lins(i)),3):Lin_Sig(lin_s(lins(i)),4);
         plot([sig(1) sig(end)],[85 85],'rx-','LineWidth',3),hold on
    end
end
if isempty(linc) == 0
    for i = 1:length(linc)
        cond = Lin_Conden(lin_c(linc(i)),4):Lin_Conden(lin_c(linc(i)),5);
        plot([cond(1) cond(end)],[80 80],'bx-','LineWidth',3)
    end
end
if isempty(poic) == 0
    for i = 1:length(poic)
        cond = Poisson_Conden(poi_c(poic(i)),4):Poisson_Conden(poi_c(poic(i)),5);
        plot([cond(1) cond(end)],[90 90],'kx-','LineWidth',3)
    end
end
if isempty(pois) == 0
    for i = 1:length(pois)
        sig = Poisson_Sig(poi_s(pois(i)),3):Poisson_Sig(poi_s(pois(i)),4);
        plot([sig(1) sig(end)],[95 95],'mx-','LineWidth',3)
    end
end
if isempty(macs) == 0
    for i = 1:length(macs)
        macwin = str2double(Macs_text(mac(macs(i)),2)):str2double(Macs_text(mac(macs(i)),3));
        plot([macwin(1) macwin(end)],[70 70],'gx-','LineWidth',3)
    end
end
tss = cell2mat(knownGene(g,6));
tes = cell2mat(knownGene(g,7));
plot(tss*ones(1,2),[0 60],'ro-',tes*ones(1,2),[0 60],'ro-','LineWidth',2)
area(win*25,Rad23b.(chr)(win)), axis([x1 x2 0 100])
set(gca,'FontSize',18,'FontWeight','bold','LineWidth',3)
ylabel('Rad23b Reads'),title(gene)
hold off
subplot(2,1,2)
area(win*25,PI.(chr)(win))
axis([x1 x2 0 100])
set(gca,'FontSize',18,'FontWeight','bold','LineWidth',2)
ylabel('PI Reads')
% subplot(3,1,3)
% area(win*25,Sox2.(chr)(win))
% axis([x1 x2 0 100])
% set(gca,'FontSize',18,'FontWeight','bold','LineWidth',2)
% ylabel('Sox2 Reads')

figure(2)
subplot(2,1,1)
win2 = floor((tss-500)/25):floor((tss+500)/25);
area(win2*25,Rad23b.(chr)(win2)),hold on
plot(win2(1)*25*ones(1,2),[0 60],'ko-',win2(end)*25*ones(1,2),[0 60],'ko-','LineWidth',2), hold off
set(gca,'FontSize',18,'FontWeight','bold','LineWidth',2),axis([x1 x2 0 100])
ylabel('Rad23b Reads'),title(gene)
subplot(2,1,2)
area(win*25,Rad23b.(chr)(win)),hold on
plot(win2(1)*25*ones(1,2),[0 60],'ko-',win2(end)*25*ones(1,2),[0 60],'ko-','LineWidth',2), hold off
set(gca,'FontSize',18,'FontWeight','bold','LineWidth',2),axis([x1 x2 0 100])
ylabel('Rad23b Reads')


