%Make all of the plots needed to visualize the PC and Zscore information

[m,n] = size(PC1);
j = 1;
for i = 1:6
    if n == 2
        x = eval(sprintf('PC%d',i));
        y = eval(sprintf('Zscore%d',i));
        figure(j)
        plot(x(:,1),x(:,2),'r*'), title(sprintf('Principle Components for Case %d',i))
        xlabel('PC 1'),ylabel('PC 2'), hold on
        figure(j+1)
        plot(y(:,1),y(:,2),'r*'), title(sprintf('Z-score values for Case %d',i))
        xlabel('TBP'),ylabel('PolII Ser5'), hold on
    else
        x = eval(sprintf('PC%d',i));
        y = eval(sprintf('Zscore%d',i));
        figure(j)
        plot3(x(:,1),x(:,2),x(:,3),'k*'), title(sprintf('Principle Components for Case %d',i))
        xlabel('PC 1'),ylabel('PC 2'),zlabel('PC 3'), hold on
        figure(j+1)
        plot3(y(:,1),y(:,3),y(:,2),'k*'), title(sprintf('Z-score values for Case %d',i))
        xlabel('TBP'),zlabel('TAF1'),ylabel('PolII Ser5'), hold on
    end
    j = j + 2;
end

