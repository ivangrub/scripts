thresh = [3 10];

for i = 1:length(thresh)
    up = VarName2./VarName3 >= thresh(i) & VarName2 >= 1;
    down = VarName3./VarName2 >= thresh(i) & VarName3 >= 1;
    figure(i)
    loglog(VarName3,VarName2,'b.',VarName3(up),VarName2(up),'r.',VarName3(down),VarName2(down),'g.'),
    hold on,plot([1e-2 1],[1 1],'k',[1 1],[1e-2 1],'k',[10^(-log10(thresh(i))) 1e3],[1 10^(3+log10(thresh(i)))],'k',[1 10^(3+log10(thresh(i)))],[10^(-log10(thresh(i))) 1e3],'k'),
    axis([1e-2 1e4 1e-2 1e4])
    xlabel('HY4')
    ylabel('HY5')
    hold off

    fid = fopen(sprintf('HY_4vs5_upgenes_%dfold.txt',thresh(i)),'w');
    for j = 1:sum(up)
        fup = find(up);
        fprintf(fid,'%s\n',char(Genes{fup(j)}));
    end
    fclose(fid);

    fid = fopen(sprintf('HY_4vs5_downgenes_%dfold.txt',thresh(i)),'w');
    for j = 1:sum(down)
        fdown = find(down);
        fprintf(fid,'%s\n',char(Genes{fdown(j)}));
    end
    fclose(fid);

    saveas(figure(i),sprintf('HY_4vs5_%dfold.pdf',thresh(i)),'pdf')
    clear fold fid
end