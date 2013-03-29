Peaks = [10138 7196;7084 3687;4895 2125;6424 1794];

[m,n] = size(Peaks);

DF = (n-1)*(m-1);

chistat = zeros(1,m);
E = zeros(1,m);

for i = 1:m
    E(i) = sum(Peaks(i,:))*sum(Peaks(:,1))/sum(sum(Peaks));
    chistat(i) = (Peaks(i,2)-E(i)).^2./E(i); 
end

1-chi2cdf(sum(chistat),DF)