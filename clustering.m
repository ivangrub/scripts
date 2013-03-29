function [idx,ctrs,sum] = clustering(Data,n)

[k,l] = size(Data);
markers = {'g*' 'r*' 'k*' 'y*' 'b*' 'm*' 'c*' 'gs' 'rs' 'ks' 'ys' 'bs' 'ms'...
    'cs' 'g^' 'r^' 'k^' 'y^' 'b^' 'm^' 'c^'};

%K-means Clustering
[idx,ctrs,sum] = kmeans(Data,n);
%subplot(1,2,1)
if l == 3
    for i = 1:n
        x = eval(sprintf('Data(idx==%d,1)',i));
        y = eval(sprintf('Data(idx==%d,2)',i));
        z = eval(sprintf('Data(idx==%d,3)',i));
        plot3(x,y,z,markers{i})
        hold on
    end
    plot3(ctrs(:,1),ctrs(:,2),ctrs(:,3),'kx','MarkerSize',12,'LineWidth',2)
    plot3(ctrs(:,1),ctrs(:,2),ctrs(:,3),'ko','MarkerSize',12,'LineWidth',2)
    xlabel('Principle Component 1'),ylabel('Principle Component 2'),zlabel('Principle Component 3')
    title('K-means Clustering')
else
    for i = 1:n
        x = eval(sprintf('Data(idx==%d,1)',i));
        y = eval(sprintf('Data(idx==%d,2)',i));
        plot(x,y,markers{i})
        hold on
    end
    plot(ctrs(:,1),ctrs(:,2),'kx','MarkerSize',12,'LineWidth',2)
    plot(ctrs(:,1),ctrs(:,2),'ko','MarkerSize',12,'LineWidth',2)
    %xlabel('TBP'),ylabel('PolII')
    %title('K-means Clustering')
end


%Fuzzy Clustering Method
%[center,U,obj_fcn] = fcm(Data,n);
%colors = {'g' 'r' 'k' 'y' 'b' 'm' 'c' 'g' 'r' 'k' 'y' 'b' 'm' 'c'...
%    'g' 'r' 'k' 'y' 'b' 'm' 'c'};
%subplot(1,2,2)  
%maxU = max(U);
%if l == 3
%    plot3(Data(:,1), Data(:,2),Data(:,3),'.')
%    for i = 1:n
%        yy = find(U(i,:) == maxU);
%        x = Data(yy,1);
%        y = Data(yy,2);
%        z = Data(yy,3);
%        if i <= 7
%            line(x,y,z,'linestyle','none','marker','*','color',colors{i})
%         elseif i > 7 && i <=14
%             line(x,y,z,'linestyle','none','marker','s','color',colors{i})
%         else
%             line(x,y,z,'linestyle','none','marker','^','color',colors{i})
%         end
%         hold on
%     end
%     plot3(center(:,1),center(:,2),center(:,3),'kx','MarkerSize',12,'LineWidth',2)
%     plot3(center(:,1),center(:,2),center(:,3),'ko','MarkerSize',12,'LineWidth',2)
%     xlabel('TBP'),ylabel('TAF1'),zlabel('PolII')
%     title('Fuzzy c-means Clustering')
% else
%    plot(Data(:,1),Data(:,2),'.')
%    for i = 1:n
%         yy = find(U(i,:) == maxU);
%         x = Data(yy,1);
%         y = Data(yy,2);
%         if i <= 7
%             line(x,y,'linestyle','none','marker','*','color',colors{i})
%         elseif i > 7 && i <=14
%             line(x,y,'linestyle','none','marker','s','color',colors{i})
%         else
%             line(x,y,'linestyle','none','marker','^','color',colors{i})
%         end
%         hold on
%     end
%     plot(center(:,1),center(:,2),'kx','MarkerSize',12,'LineWidth',2)
%     plot(center(:,1),center(:,2),'ko','MarkerSize',12,'LineWidth',2)
%     xlabel('TBP'),ylabel('PolII')
%     title('Fuzzy c-means Clustering')
% end

end