%Make the average reads plot for the TSS

win = 5000;
avgread = zeros((win+read_length-1)*2+1,length(headers)+1);
for i = 1:length(headers)
    x = eval(headers{i});
    avg = bp2TSS(x,knownGene,read_length,win);
    avgread(:,[1 i+1]) = avg;
end

save(sprintf('%s_SingleBP_10000bp_%s.mat',perID{1},date),'ReadSum','headers','avgread')
