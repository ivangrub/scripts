%TAF7 and Rad23b Coverage

filename = 'OverlapPeaksFor_Rad23b_relativeto_TAF7.txt';

fid = fopen(filename);
data = textscan(fid,'%s%d%d%s%d%d%d%d%d%d','HeaderLines',1);
fclose(fid);clear fid filename

cover = ~strcmp(data{1,4},'NaN');
Overlapped = sum(cover);

Peaks = find(cover);

win = zeros(4,sum(cover));

win(1,:) = (data{1,5}(cover) >= data{1,2}(cover) & data{1,6}(cover) <= data{1,3}(cover));       %IP window fits inside of OS window
win(2,:) = (data{1,5}(cover) < data{1,2}(cover) & data{1,6}(cover) >= data{1,2}(cover));       %IP window straddles start of OS window
win(3,:) = (data{1,5}(cover) <= data{1,3}(cover) & data{1,6}(cover) > data{1,3}(cover));       %IP window staddles the end of the window
win(4,:) = (data{1,5}(cover) <= data{1,2}(cover) & data{1,6}(cover) >= data{1,3}(cover));       %OS window fits inside of the IP window

Coverage = zeros(1,sum(cover));

for i = 1:length(Peaks)
    if win(1,i)
        Coverage(i) = (double(data{1,6}(Peaks(i))) - double(data{1,5}(Peaks(i))))/(double(data{1,3}(Peaks(i))) - double(data{1,2}(Peaks(i))));
        continue
    end
    if win(2,i)
        Coverage(i) = (double(data{1,6}(Peaks(i))) - double(data{1,2}(Peaks(i))))/(double(data{1,3}(Peaks(i))) - double(data{1,2}(Peaks(i))));
    end
    if win(3,i)
        Coverage(i) = (double(data{1,3}(Peaks(i))) - double(data{1,5}(Peaks(i))))/(double(data{1,3}(Peaks(i))) - double(data{1,2}(Peaks(i))));
    end
    if win(4,i)
        Coverage(i) = 1;
        continue
    end
end

[freq,xout]=hist(Coverage,100);
plot(xout,freq/length(Coverage),'r')