%A plot of the average number of reads per window as a function of distanceance
%from the TSS

function [ys ye] = AverageReadsPlot(distance,Reads,headers)

[n,~] = size(Reads);
x = (-5:.25:5)*1000;
ys = zeros(length(x)-1,n);
ye = zeros(length(x)-1,n);

for i = 1:length(x)-1
    a = (distance(:,1) >= x(i) & distance(:,1) < x(i+1)) | (distance(:,2) >= x(i) & distance(:,2) < x(i+1));
    b = (distance(:,3) >= x(i) & distance(:,3) < x(i+1)) | (distance(:,4) >= x(i) & distance(:,4) < x(i+1));
    ys(i,:) = mean(Reads(a,:),1);
    ye(i,:) = mean(Reads(b,:),1);
end

plot(x(2:end),ys)
xlabel('Distance from TSS/TES (bp)')
ylabel('Average Reads per Window')
legend(headers),hold on
plot(x(2:end),ye), hold off


clear i x ys ye n a b c d              