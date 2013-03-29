subplot(1,3,1)
bar(ESCCoeff1(:,1)), xlabel('1st Principal Component')
subplot(1,3,2)
bar(ESCCoeff1(:,2)), xlabel('2nd Principal Component')
subplot(1,3,3)
bar(ESCCoeff1(:,3)), xlabel('3rd Principal Component')


tbp = zscore(chip(:,1));
%taf = zscore(chip(:,2));
pol = zscore(chip(:,2));

tbplog = zscore(loggedchip(:,1));
%taflog = zscore(loggedchip(:,2));
pollog = zscore(loggedchip(:,2));

figure(1)
subplot(3,1,1)
plot(tbp,taf,'r*'), hold on
plot(tbp(49410:end),taf(49410:end),'b*'),hold off
xlabel('TBP Z-score'),ylabel('TAF1 Z-score')
legend('Data','tRNA')
title('ESC Z-scores of Normalized Data for Linear Scale')
axis square
subplot(3,1,2)
plot(tbp,pol,'r*'), hold on
plot(tbp(49410:end),pol(49410:end),'b*'),hold off
xlabel('TBP Z-score'),ylabel('Pol2 Z-score')
legend('Data','tRNA')
axis square
subplot(3,1,3)
plot(taf,pol,'r*')
xlabel('TAF1 Z-score'),ylabel('Pol2 Z-score')
axis square

figure(2)
subplot(3,1,1)
plot(tbplog,taflog,'r*'), hold on
plot(tbplog(49410:end),taflog(49410:end),'b*'),hold off
xlabel('TBP Z-score'),ylabel('TAF1 Z-score')
legend('Data','tRNA')
title('ESC Z-scores of Normalized Data for Log Scale')
subplot(3,1,2)
plot(tbplog,pollog,'r*'), hold on
plot(tbplog(49410:end),pollog(49410:end),'b*'),hold off
xlabel('TBP Z-score'),ylabel('Pol2 Z-score')
legend('Data','tRNA')
axis square
subplot(3,1,3)
plot(taflog,pollog,'r*')
xlabel('TAF1 Z-score'),ylabel('Pol2 Z-score')
axis square

tbpreads = sum(CellType(:,1));
tafreads = sum(CellType(:,3));
polreads = sum(CellType(:,5));

normed = normalize(CellType(:,1:2:9));
x = zscore(normed(:,1),1);
y = zscore(cpgchip(:,1),1);
j = find(x >= 3);
k = find(y >= 3);
subplot(2,2,1)
histfit(normed(:,1),50)
xlabel('Number of Reads for TBP'),ylabel('Number of Occurances')
title('Linear Scale of Enrichment')
subplot(2,2,3)
histfit(cpgchip(:,1),50)
xlabel('Number of Reads for TBP'),ylabel('Number of Occurances')
title('Log Scale of Enrichment')
subplot(2,2,2)
hist(x,50)
xlabel('Z-score of TBP Enrichment'),ylabel('Number of Occurances')
title('Linear Scale of Enrichment')
subplot(2,2,4)
hist(y,50)
xlabel('Z-score of TBP Enrichment'),ylabel('Number of Occurances')
title('Log Scale of Enrichment')
m = x >= 3;
n = y >=3;
sum(m)
sum(n)


[PCtbppol Zscoretbppol coefftbppol latent cumvar] = PCA(chip(:,[1 3]));

subplot(1,2,1)
bar(coefftbppol(:,1)), xlabel('1st Principal Component')
subplot(1,2,2)
bar(coefftbppol(:,2)), xlabel('2nd Principal Component')
title('PCA of ESC for TBP and Pol2')

figure(1)
subplot(1,2,1)
plot(tbp,pol,'r*'), hold on
plot(tbp(49410:end),pol(49410:end),'b*'),hold off
xlabel('TBP Z-score'),ylabel('Pol2 Z-score')
legend('Data','tRNA')
title('Motor Neuron Z-scores on a Linear Scale')
subplot(1,2,2)
plot(tbplog,pollog,'r*'), hold on
plot(tbplog(49410:end),pollog(49410:end),'b*'),hold off
xlabel('TBP Z-score'),ylabel('Pol2 Z-score')
legend('Data','tRNA')
title('Motor Neuron Z-scores on a Log Scale')

subplot(1,3,1)
plot(ESCZscore1(49410:end,1),ESCZscore1(49410:end,2),'r*'),hold on
plot(ESCZscore1(28645,1),ESCZscore1(28645,2),'k*',ESCZscore1(28647,1),ESCZscore1(28647,2),'k*',...
    ESCZscore1(28650,1),ESCZscore1(28650,2),'k*',ESCZscore1(28652,1),ESCZscore1(28652,2),'k*',...
    ESCZscore1(28653,1),ESCZscore1(28653,2),'k*')
plot(ESCZscore1(5482,1),ESCZscore1(5482,2),'c*',ESCZscore1(5985,1),ESCZscore1(5985,2),'c*',...
    ESCZscore1(6057,1),ESCZscore1(6057,2),'c*',ESCZscore1(6341,1),ESCZscore1(6341,2),'c*',...
    ESCZscore1(6865,1),ESCZscore1(6865,2),'c*')
plot(ESCZscore1(37541,1),ESCZscore1(37541,3),'g*')
plot(ESCZscore1(19048,1),ESCZscore1(19048,3),'b*')
plot(ESCZscore1(27585,1),ESCZscore1(27585,3),'m*')
xlabel('TBP Z-score'),ylabel('TAF1 Z-score'),axis square, grid on
subplot(1,3,2)
plot(ESCZscore1(49410:end,1),ESCZscore1(49410:end,3),'r*'),hold on
plot(ESCZscore1(28645,1),ESCZscore1(28645,3),'k*',ESCZscore1(28647,1),ESCZscore1(28647,3),'k*',...
    ESCZscore1(28650,1),ESCZscore1(28650,3),'k*',ESCZscore1(28652,1),ESCZscore1(28652,3),'k*',...
    ESCZscore1(28653,1),ESCZscore1(28653,3),'k*')
plot(ESCZscore1(5482,1),ESCZscore1(5482,3),'c*',ESCZscore1(5985,1),ESCZscore1(5985,3),'c*',...
    ESCZscore1(6057,1),ESCZscore1(6057,3),'c*',ESCZscore1(6341,1),ESCZscore1(6341,3),'c*',...
    ESCZscore1(6865,1),ESCZscore1(6865,3),'c*')
plot(ESCZscore1(37541,1),ESCZscore1(37541,3),'g*')
plot(ESCZscore1(19048,1),ESCZscore1(19048,3),'b*')
plot(ESCZscore1(27585,1),ESCZscore1(27585,3),'m*')
xlabel('TBP Z-score'),ylabel('TAF1 Z-score'),axis square, grid on
subplot(1,3,3)
plot(ESCZscore1(49410:end,2),ESCZscore1(49410:end,3),'r*'),hold on
plot(ESCZscore1(28645,2),ESCZscore1(28645,3),'k*',ESCZscore1(28647,2),ESCZscore1(28647,3),'k*',...
    ESCZscore1(28650,2),ESCZscore1(28650,3),'k*',ESCZscore1(28652,2),ESCZscore1(28652,3),'k*',...
    ESCZscore1(28653,2),ESCZscore1(28653,3),'k*')
plot(ESCZscore1(5482,2),ESCZscore1(5482,3),'c*',ESCZscore1(5985,2),ESCZscore1(5985,3),'c*',...
    ESCZscore1(6057,2),ESCZscore1(6057,3),'c*',ESCZscore1(6341,2),ESCZscore1(6341,3),'c*',...
    ESCZscore1(6865,2),ESCZscore1(6865,3),'c*')
plot(ESCZscore1(37541,2),ESCZscore1(37541,3),'g*')
plot(ESCZscore1(19048,2),ESCZscore1(19048,3),'b*')
plot(ESCZscore1(27585,2),ESCZscore1(27585,3),'m*')
xlabel('TBP Z-score'),ylabel('TAF1 Z-score')axis square, grid on

biplot(ESCCoeff1(:,1:3),'Scores',ESCPCTotal1(:,1:3))

subplot(1,3,1)
plot(CellType(:,1),CellType(:,3),'*')
xlabel('Normalized TBP Enrichment'),ylabel('Normalized TAF1 Enrichment')
axis square
subplot(1,3,2)
plot(CellType(:,1),CellType(:,5),'*')
xlabel('Normalized TBP Enrichment'),ylabel('Normalized Pol2 Enrichment')
axis square
subplot(1,3,3)
plot(CellType(:,3),CellType(:,5),'*')
xlabel('Normalized TAF1 Enrichment'),ylabel('Normalized Pol2 Enrichment')
axis square 
subplot(2,1,2)
plot(tafz,ntafz,'*'),xlabel('TAF1 Enrichment'), ylabel('Normalized TAF1 Enrichment')

subplot(2,2,1)
plot(ntaf,taf,'*'),xlabel('Normalized TAF1 Enrichment'),ylabel('TAF1 Enrichment')
subplot(2,2,2)

