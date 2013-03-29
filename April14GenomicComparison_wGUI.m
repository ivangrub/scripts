fh=figure;
pbh1 = uicontrol(fh,'Style','pushbutton','String','Internal',...
                'Position',[10 300 60 40]);
pbh2 = uicontrol(fh,'Style','pushbutton','String','Pre-Immune',...
                'Position',[10 250 60 40]);
pbh3 = uicontrol(fh,'Style','pushbutton','String','Macs',...
                'Position',[10 200 60 40]);
pbh4 = uicontrol(fh,'Style','pushbutton','String','Go',...
                'Position',[400 350 60 40]);
gene = uicontrol(fh,'Style','edit',...
                'String','Enter Gene Name',...
                'Position',[225 400 130 20]);

if pbh4 == 1
g = find(strcmp(gene,Lin_Conden_text(:,3)),1,'first');
chr = Lin_Conden_text(g,5);
x1 = str2double(Lin_Conden_text(g,7))-8000;
x2 = str2double(Lin_Conden_text(g,8))+8000;


win = floor(x1/25):floor(x2/25);


if pbh2 == 1
    lin_c = find(strcmp(chr,Lin_Conden_text(:,5)));
    linc = find(Lin_Conden(lin_c,4) > x1 & Lin_Conden(lin_c,5) < x2);
    lin_s = find(strcmp(chr,Lin_Sig_text(:,5)));
    lins = find(Lin_Sig(lin_s,3) > x1 & Lin_Sig(lin_s,4) < x2);
    if isempty(lins) == 0
        for i = 1:length(lins)
             sig = Lin_Sig(lin_s(lins(i)),3):Lin_Sig(lin_s(lins(i)),4);
             plot([sig(1) sig(end)],[85 85],'rx-','LineWidth',3), hold on
        end
    end

    if isempty(linc) == 0
        for i = 1:length(linc)
            cond = Lin_Conden(lin_c(linc(i)),4):Lin_Conden(lin_c(linc(i)),5);
            plot([cond(1) cond(end)],[80 80],'bx-','LineWidth',3)
        end
    end
end
if pbh1 == 1
    poi_c = find(strcmp(chr,Poisson_Conden_text(:,5)));
    poic = find(Poisson_Conden(poi_c,4) > x1 & Poisson_Conden(poi_c,5) < x2);
    poi_s = find(strcmp(chr,Poisson_Sig_text(:,5)));
    pois = find(Poisson_Sig(poi_s,3) > x1 & Poisson_Sig(poi_s,4) < x2);
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
end
if pbh3 == 1
    mac = find(strcmp(chr,Macs_text(:,1)));
    macs = find(str2double(Macs_text(mac,2)) > x1 & str2double(Macs_text(mac,3)) < x2);
    if isempty(macs) == 0
        for i = 1:length(macs)
            macwin = str2double(Macs_text(mac(macs(i)),2)):str2double(Macs_text(mac(macs(i)),3));
            plot([macwin(1) macwin(end)],[70 70],'gx-','LineWidth',3)
        end
    end
end
plot(str2double(Lin_Conden_text(g,7))*ones(1,2),[0 60],'ro-',str2double(Lin_Conden_text(g,8))*ones(1,2),[0 60],'ro-','LineWidth',2)
area(win*25,Rad23b.(chr{1})(win)), axis([x1 x2 0 100])
set(gca,'FontSize',18,'FontWeight','bold','LineWidth',2)
ylabel('Rad23b Reads')
hold off
end
% subplot(3,1,2)
% area(win*25,PI.(chr{1})(win))
% axis([x1 x2 0 100])
% set(gca,'FontSize',18,'FontWeight','bold','LineWidth',2)
% ylabel('PI Reads')
% 
% subplot(3,1,3)
% win2 = floor((str2double(Lin_Conden_text(g,7))-500)/25):floor((str2double(Lin_Conden_text(g,7))+500)/25);
% area(win*25,Rad23b.(chr{1})(win)),hold on
% plot(win2(1)*25*ones(1,2),[0 60],'ko-',win2(end)*25*ones(1,2),[0 60],'ko-','LineWidth',2), hold off
% set(gca,'FontSize',18,'FontWeight','bold','LineWidth',2),axis([x1 x2 0 100])
% ylabel('Rad23b Reads')



