for i = 1:8
    [idx x sumd] = kmeans(ESCPCTotal1,i);
    l(i) = sum(sumd);
end
figure(1)
plot(l), xlabel('Number of Clusters'),ylabel('Total Distance to the Centroid')

[idx x sumd] = kmeans(ESCPCTotal1,4);
   
figure(1)
subplot(4,1,1)
bar(mean(ESCZscore1(idx == 1,1:3),1)'), hold on
errorbar(mean(ESCZscore1(idx == 1,1:3),1)',std(ESCZscore1(idx == 1,1:3),0,1)','xk')
ylabel('Cluster 1')
title('Z-scores for the Log Scaled Enrichment Populations for each Cluster')
hold off
subplot(4,1,2)
bar(mean(ESCZscore1(idx == 2,1:3),1)), hold on
errorbar(mean(ESCZscore1(idx == 2,1:3),1),std(ESCZscore1(idx == 2,1:3),0,1),'xk')
ylabel('Cluster 2')
hold off
subplot(4,1,3)
bar(mean(ESCZscore1(idx == 3,1:3),1)), hold on
errorbar(mean(ESCZscore1(idx == 3,1:3),1),std(ESCZscore1(idx == 3,1:3),0,1),'xk')
ylabel('Cluster 3')
hold off
subplot(4,1,4)
bar(mean(ESCZscore1(idx == 4,1:3),1)), hold on
errorbar(mean(ESCZscore1(idx == 4,1:3),1),std(ESCZscore1(idx == 4,1:3),0,1),'xk')
ylabel('Cluster 4'), xlabel('Transcription Factors TBP, TAF1 and Pol2 Respectively')
hold off
% subplot(5,1,5)
% bar(mean(ESCZscore1(idx == 5,1:3),1)), hold on
% errorbar(mean(ESCZscore1(idx == 5,1:3),1),std(ESCZscore1(idx == 5,1:3),0,1),'xk')
% hold off


figure(2)
plot3(ESCZscore1(idx==1,1),ESCZscore1(idx==1,2),ESCZscore1(idx==1,3),'r*'), hold on
plot3(ESCZscore1(idx==2,1),ESCZscore1(idx==2,2),ESCZscore1(idx==2,3),'g*'),
plot3(ESCZscore1(idx==3,1),ESCZscore1(idx==3,2),ESCZscore1(idx==3,3),'k*'),
plot3(ESCZscore1(idx==4,1),ESCZscore1(idx==4,2),ESCZscore1(idx==4,3),'b*'),
%plot3(ESCZscore1(idx==5,1),ESCZscore1(idx==5,2),ESCZscore1(idx==5,3),'y*'),
xlabel('TBP Z-score'),ylabel('TAF1 Z-score'),zlabel('Pol2 Z-score')
legend('Cluster 1','Cluster 2' ,'Cluster 3', 'Cluster 4')
title('Clusters for Log Scale ESC Data')

