Peaks = [10138 2048;7084 1899;4895 867;6424 535];

combos = [1 2;1 3;1 4;2 3;2 4;3 4];

[m,~] = size(combos);

z = zeros(1,m);
for i = 1:m
    pstat = (Peaks(combos(i,1),2) + Peaks(combos(i,2),2))/(Peaks(combos(i,1),1) + Peaks(combos(i,2),1));
    p1 = Peaks(combos(i,1),2)/Peaks(combos(i,1),1);
    p2 = Peaks(combos(i,2),2)/Peaks(combos(i,2),1);
    n1 = Peaks(combos(i,1),1);
    n2 = Peaks(combos(i,2),1);
    z(i) = (p1-p2)/sqrt(pstat*(1-pstat)*(1/n1+1/n2));
end

double(1-normcdf(z,0,1))