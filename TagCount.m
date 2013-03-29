%Count the full read count

ReadSum = zeros(1,length(headers));
for i = 1:length(headers)
    x = eval(headers{i});
    chr = fieldnames(x.win);
    for j = 1:length(chr)
        ReadSum(i) = ReadSum(i) + sum(x.win.(chr{j})(2,:));
    end
end