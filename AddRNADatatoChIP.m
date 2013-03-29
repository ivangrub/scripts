clear,clc

load FPKM_ESC.mat
headers = {'Rad23b_Original' 'Oct4' 'Sox2'};
comment = 'Lin';

for j = 1:length(headers)
    x = importdata(sprintf('%s_%s_Conden.txt',headers{j},comment));
    x.data(:,9:10) = zeros(length(x.data),2);
    for i = 2:length(x.textdata)
        y = find(strcmp(x.textdata(i,1),FPKM_ESC(:,1)));
        if ~isempty(y) && ~isempty(FPKM_ESC{y,4})
            x.data(i-1,9) = FPKM_ESC{y,4}(1);
            x.data(i-1,10) = FPKM_ESC{y,5}(1);
        else
            x.data(i-1,9:10) = [NaN NaN];
        end
    end
    assignin('base',sprintf('%s_Conden_wRNA',headers{j}),x.data);
    assignin('base',sprintf('%s_Conden_wRNA_text',headers{j}),x.textdata);
end

