plot(ESCZscore1(49410:end,1),ESCZscore1(49410:end,2),'r*'),hold on
plot(ESCZscore1(28645,1),ESCZscore1(28645,2),'k*',ESCZscore1(28647,1),ESCZscore1(28647,2),'k*',...
ESCZscore1(28650,1),ESCZscore1(28650,2),'k*',ESCZscore1(28652,1),ESCZscore1(28652,2),'k*',...
ESCZscore1(28653,1),ESCZscore1(28653,2),'k*')
plot(ESCZscore1(5482,1),ESCZscore1(5482,2),'c*',ESCZscore1(5985,1),ESCZscore1(5985,2),'c*',...
ESCZscore1(6057,1),ESCZscore1(6057,2),'c*',ESCZscore1(6341,1),ESCZscore1(6341,2),'c*',...
ESCZscore1(6865,1),ESCZscore1(6865,2),'c*')
plot(ESCZscore1(37541,1),ESCZscore1(37541,2),'g*')
plot(ESCZscore1(19048,1),ESCZscore1(19048,2),'b*')
plot(ESCZscore1(27585,1),ESCZscore1(27585,2),'y*')
xlabel('TBP Z-score'),ylabel('TAF1 Z-score'),axis square, grid on

subplot(1,2,1)
plot(ESCZscore1(:,1),ESCZscore1(:,2),'b*',ESCZscore1(49410:end,1),ESCZscore1(49410:end,2),'r*',...
    ESCZscore1(19048,1),ESCZscore1(19048,2),'g*',ESCZscore1(27585,1),ESCZscore1(27585,2),'k*')
xlabel('TBP Z-score'),ylabel('TAF1 Z-score'),axis square
subplot(1,2,2)
plot(ESCZscore1(:,3),ESCZscore1(:,2),'b*',ESCZscore1(49410:end,3),ESCZscore1(49410:end,2),'r*',...
    ESCZscore1(19048,3),ESCZscore1(19048,2),'g*',ESCZscore1(27585,3),ESCZscore1(27585,2),'k*')
xlabel('Pol2 Z-score'),ylabel('TAF1 Z-score'),legend('Remaining Data','tRNA','Oct4','Sox2'),axis square





plot3(ESCZscore1(:,1),ESCZscore1(:,2),ESCZscore1(:,3),'b*',ESCZscore1(c5,1),ESCZscore1(c5,2),ESCZscore1(c5,3),'r*',...
    ESCZscore1(c7,1),ESCZscore1(c7,2),ESCZscore1(c7,3),'k*',...
    ESCZscore1(c4,1),ESCZscore1(c4,2),ESCZscore1(c4,3),'g*')
xlabel('TBP Z-score'),ylabel('TAF1 Z-score'),zlabel('Pol2 Z-score'),axis square, grid on
legend('All Data','Only TBP Cluster','TBP and TAF1 Cluster','TAF1 and Pol2 Cluster')

Genescores = dot(ESCPCTotal1,[-1 -1 1]);