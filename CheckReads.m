%Subdivisions to accept as a positive hit. Ranges that will be inspected
%will be as follows: [3 4] [4 8] [8 12]
%Anything less than z = 1.5 will be considered to be empty.

%lim = [3 4 6 12];
%gene = {'uniquegene' 'cpggene' 'restgene'};
k = 1;

x = ESCZscore1;
es1 = (x(:,1) >= 3 & x(:,2) <= 1.5 & x(:,3) <= 1.5);
es2 = (x(:,2) >= 3 & x(:,1) <= 1.5 & x(:,3) <= 1.5);
es3 = (x(:,3) >= 3 & x(:,2) <= 1.5 & x(:,1) <= 1.5);


if isempty(es1) ~= 1
    I = find(es1);
    for z = 1:5:length(I)
        j = I(z);
        if isnan(str2double(knownGene{j,4}))
            TBP = eval(sprintf('tbpchip.chr%s',knownGene{j,4}));
            TAF = eval(sprintf('tafchip.chr%s',knownGene{j,4}));
            POL2 = eval(sprintf('polchip.chr%s',knownGene{j,4}));
            PI = eval(sprintf('pichip.chr%s',knownGene{j,4}));
        else
            TBP = eval(sprintf('tbpchip.chr%d',str2double(knownGene{j,4})));
            TAF = eval(sprintf('tafchip.chr%d',str2double(knownGene{j,4})));
            POL2 = eval(sprintf('polchip.chr%d',str2double(knownGene{j,4})));
            PI = eval(sprintf('pichip.chr%d',str2double(knownGene{j,4})));
        end
        figure(1)
        plot(TBP(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'r'),hold on
        plot(TAF(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'k')
        plot(POL2(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'b')
        plot(PI(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'g')
        ylabel('Number of Reads')
        xlabel('1000 Base Pair Window')
        legend('TBP','TAF1','PolII','PI')
        title(sprintf('%s-%s  is TBP specific',knownGene{j,2},knownGene{j,3}))
        hold off
        saveas(figure(1),sprintf('checkreads%d',k),'png')
        k = k + 1;
    end
end
if isempty(es2) ~= 1
    I = find(es2);
    for z = 1:5:length(I)
        j = I(z);
        if isnan(str2double(knownGene{j,4}))
            TBP = eval(sprintf('tbpchip.chr%s',knownGene{j,4}));
            TAF = eval(sprintf('tafchip.chr%s',knownGene{j,4}));
            POL2 = eval(sprintf('polchip.chr%s',knownGene{j,4}));
            PI = eval(sprintf('pichip.chr%s',knownGene{j,4}));
        else
            TBP = eval(sprintf('tbpchip.chr%d',str2double(knownGene{j,4})));
            TAF = eval(sprintf('tafchip.chr%d',str2double(knownGene{j,4})));
            POL2 = eval(sprintf('polchip.chr%d',str2double(knownGene{j,4})));
            PI = eval(sprintf('pichip.chr%d',str2double(knownGene{j,4})));
        end
        figure(1)
        plot(TBP(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'r'),hold on
        plot(TAF(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'k')
        plot(POL2(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'b')
        plot(PI(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'g')
        ylabel('Number of Reads')
        xlabel('1000 Base Pair Window')
        legend('TBP','TAF1','PolII','PI')
        title(sprintf('%s-%s is TAF specific',knownGene{j,2},knownGene{j,3}))
        hold off
        saveas(figure(1),sprintf('checkreads%d',k),'png')
        k = k + 1;
    end
end
if isempty(es3) ~= 1
   I = find(es3);
    for z = 1:5:length(I)
        j = I(z);
        if isnan(str2double(knownGene{j,4}))
            TBP = eval(sprintf('tbpchip.chr%s',knownGene{j,4}));
            TAF = eval(sprintf('tafchip.chr%s',knownGene{j,4}));
            POL2 = eval(sprintf('polchip.chr%s',knownGene{j,4}));
            PI = eval(sprintf('pichip.chr%s',knownGene{j,4}));
        else
            TBP = eval(sprintf('tbpchip.chr%d',str2double(knownGene{j,4})));
            TAF = eval(sprintf('tafchip.chr%d',str2double(knownGene{j,4})));
            POL2 = eval(sprintf('polchip.chr%d',str2double(knownGene{j,4})));
            PI = eval(sprintf('pichip.chr%d',str2double(knownGene{j,4})));
        end
        figure(1)
        plot(TBP(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'r'),hold on
        plot(TAF(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'k')
        plot(POL2(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'b')
        plot(PI(floor(knownGene{j,6}/25-20):floor(knownGene{j,6}/25+20)),'g')
        ylabel('Number of Reads')
        xlabel('1000 Base Pair Window')
        legend('TBP','TAF1','PolII','PI')
        title(sprintf('%s-%s is Pol2 specific',knownGene{j,2},knownGene{j,3}))
        hold off
        saveas(figure(1),sprintf('checkreads%d',k),'png')
        k = k + 1;
    end
end

