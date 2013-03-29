coor = [91439302 91465882;55362914 55405109;70158904 70163683];
genes = {'XPC' 'Rad23b' 'Cetn2'};
chr = {'chr6' 'chr4' 'chrX'};
plusminus = [-100e3 100e3];

[m,~] = size(coor);

l = 1;
for j = 1:length(headers)
    x = eval(headers{j});
    for i = 1:m
    	c = round((coor(i,1)+plusminus(1))/window):round((coor(i,2)+plusminus(2))/window);
    	tss = coor(i,1);
    	tes = coor(i,2);
    	top = max(x.win.(chr{i})(2,c));
    	figure(l)
    	
    	plot(tss*ones(1,2),[0 top],'ro-',tes*ones(1,2),[0 top],'ro-','LineWidth',2),hold on
    	area(x.win.(chr{i})(1,c),x.win.(chr{i})(2,c))
        if top == 0
            axis([x.win.(chr{i})(1,c(1)) x.win.(chr{i})(1,c(end)) 0 1]);
        else
            axis([x.win.(chr{i})(1,c(1)) x.win.(chr{i})(1,c(end)) 0 1.2*top]);
        end
        title(sprintf('%s on %s',headers{j},genes{i})), hold off
        l = l + 1;
        saveas(gcf,sprintf('Yick_Marson_CellPaper_%s_on_%s_woSPP.eps',headers{j},genes{i}));
    end
end