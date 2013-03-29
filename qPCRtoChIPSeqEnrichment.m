coor = [143668352 143668484];

chr = {'chr5'};

[m,~] = size(coor);

A = zeros(length(headers),m);
B = zeros(length(headers),m);

for j = 1:length(headers)
    x = eval(headers{j});
    for i = 1:m
        A(j,i) = sum(x.bp.(chr{i}) >= coor(i,1) & x.bp.(chr{i}) <= coor(i,2));
    end
    B(j,:) = A(j,:)./ReadSum(j);
end