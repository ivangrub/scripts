function Chr = chr_save(win,data1,mm9length)

m = length(fieldnames(data1.chip));
yy = strrep(fieldnames(data1.chip),'chr','');
chr = cellstr(yy);
k = 1;
for i = 1:m;
    x = floor(mm9length(i)/25);
    tss1 = 1:win:x-win+1;
    y = length(tss1);
    B(k:y+k-1) = chr(i);
    k = k+y;
end
Chr = B';

clear m yy chr i x tss1 y k B

end