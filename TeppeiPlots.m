wTES = 0;
axisoption = 1;
headersind = [1 2 3 4];


Genes = {'Gapdh' 'Pou5f1' 'Cnn2' 'Actb' 'Sox2'};

top = zeros(1,length(headersind));
for i = 1:length(Genes)
    for j = 1:length(headersind)
        x = eval(headers{headersind(j)});
        c = find(strcmp(Genes{i},knownGene(:,2)),1,'first');
        chr = knownGene(c,4);
        if isempty(chr)
            error('%s is an incorrect gene abbreviation',Genes{i})
        end
        if cell2mat(knownGene(c,5)) == 1
            TSS = cell2mat(knownGene(c,6));
            TES = cell2mat(knownGene(c,7));
            strand = '+';
        else
            TSS = cell2mat(knownGene(c,7));
            TES = cell2mat(knownGene(c,6));
            strand = '-';
        end
        
        %Set indices and find the max peak for axis purposes
        ind = round(TSS/150);
        
        
        %Make the figure
        
        %w/o TES
        if wTES == 0
            top(j) = max(x.win.(chr{1})(2,ind-10:ind+10));
            figure(i)
            set(gcf,'Name',sprintf('%s: %s strand',Genes{i},strand))
            subplot(length(headersind),1,j)
            area(x.win.(chr{1})(1,ind-10:ind+10),x.win.(chr{1})(2,ind-10:ind+10)),hold on
            plot([TSS;TSS],[0;top(j)],'g-','LineWidth',3),hold off
            set(gca,'LineWidth',2,'FontSize',13,'FontWeight','bold')
            ylabel(test{headersind(j)})
        else
            figure(i)
            set(gcf,'Name',sprintf('%s: %s strand',Genes{i},strand))
            subplot(length(headersind),1,j)
            if strand == '+'
                top(j) = max(x.win.(chr{1})(2,ind-5:floor(TES/150)+5));
                area(x.win.(chr{1})(1,ind-5:round(TES/150)+5),x.win.(chr{1})(2,ind-5:round(TES/150)+5)),hold on
            else
                top(j) = max(x.win.(chr{1})(2,floor(TES/150)-5:ind+5));
                area(x.win.(chr{1})(1,round(TES/150)-5:ind+5),x.win.(chr{1})(2,round(TES/150)-5:ind+5)),hold on
            end
            plot([TSS;TSS],[0;top(j)],'g-',[TES;TES],[0;top(j)],'r-','LineWidth',3),hold off
            set(gca,'LineWidth',2,'FontSize',13,'FontWeight','bold')
            ylabel(test{headersind(j)})
        end
    end
    subplot(length(headersind),1,1)
    title(sprintf('%s: %s strand on %s',Genes{i},strand,chr{1}))
    
    %Re-adjust for the axis
    if axisoption == 1
        if wTES == 0
            [~,j] = max(top);
            for k = 1:length(headersind)
                subplot(length(headersind),1,k)
                axis([x.win.(chr{1})(1,ind-10) x.win.(chr{1})(1,ind+10) 0 1.2*top(j)]),hold off
            end
        else
            if strand == '+'
                [~,j] = max(top);
                for k = 1:length(headersind)
                    subplot(length(headersind),1,k)
                    axis([x.win.(chr{1})(1,ind-5) x.win.(chr{1})(1,round(TES/150)+5) 0 1.2*top(j)]),hold off
                end
            else
                [~,j] = max(top);
                for k = 1:length(headersind)
                    subplot(length(headersind),1,k)
                    axis([x.win.(chr{1})(1,round(TES/150)-5) x.win.(chr{1})(1,ind+5) 0 1.2*top(j)]),hold off
                end
            end
        end
    end
    
    if wTES == 1
        saveas(gcf,sprintf('TeppeiFigures/%s_wTES_%s_axis%d.eps',Genes{i},num2str(headersind),axisoption),'psc2');
    else
        saveas(gcf,sprintf('TeppeiFigures/%s_onlyTSS_%s_axis%d.eps',Genes{i},num2str(headersind),axisoption),'psc2');
    end
    
end

