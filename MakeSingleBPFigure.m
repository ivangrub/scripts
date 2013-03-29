%Take the SingleBP mat file and make a figure

perID = 'HY';

win = 1000;
headersind = [1 2 5 7 9 11];
test = strrep(headers,'_',' ');
y = abs(avgread(:,1)) <= win;

if win > 5000
    error('The data is limited to a 5,000 bp window')
end

x = zeros(sum(y),length(headersind));
for i = 1:length(headersind)
    x(:,i) = avgread(y,headersind(i)+1)./ReadSum(headersind(i))*100;
end

plot(avgread(y,1),x,'LineWidth',2)
legend(test(headersind))
xlabel('Distance from TSS')
ylabel('Average % of Total Read Coverage per BP')

saveas(gcf,sprintf('%s_AverageCoverage_%dbp_of_TSS_%s.pdf',perID,win,num2str(headersind)));