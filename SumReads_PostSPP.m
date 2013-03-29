%Calculate the read coverage across the whole genome for single bp and
%condensed windows

clear,clc

%% Add options and identify people
person = {'Smith' 'Smith'};
headers = {'PI' 'ELL'};
test = strrep(headers,'_',' ');
window = 250;
read_length = 50;
dna_len = 250;

perID = {'Smith'};

data = importdata('mm9.chr.length');
len = data.data;
chr = data.textdata;
clear data

%% Initiate matrices
for i = 1:length(headers)
    for j = 1:length(chr)
        bp.(chr{j}) = zeros(1,round(len(j)/window));
        left = 1:window:len(j)-window;
        right = zeros(1,round(len(j)/window));
        if length(left) ~= length(right)
            if length(left) > length(right)
                left = left(1:length(right));
            else
                left(length(left)+1) = left(end) + window;
            end
        end
        win.(chr{j}) = [left;right];
    end
    Y = struct('bp',bp,'win',win);
    assignin('base',headers{i},Y);
    clear bp win Y  
end

clear len chr
%% Scan through the data
for i = 1:length(headers)
    fprintf(sprintf('Working on %s\n',headers{i}))
    fid = fopen(sprintf('~/ChIP-Data/%s/%s_%s_wSPP.tags.txt',person{1},person{2},headers{i}),'r');
    x = eval(headers{i});
    y = textscan(fid,'%d%d%d%d%d','HeaderLines',1);
    fclose(fid);clear fid
    B = [y{1,1} y{1,2} y{1,3} y{1,4} y{1,5}]';
    
    A = reshape(B,5*length(y{1,1}),1)';
    
    fid = fopen(sprintf('~/ChIP-Data/%s/%s_%s_SPP_chr.txt',person{1},person{2},headers{i}),'r');
    XX = textscan(fid,'%s');
    fclose(fid);clear fid
    CHR = XX{1,1};
    
    clear y z f XX B
    
    if isempty(A)
        error('There is an error with loading the SPP Tags File')
    end
    
    j = 1;
    len = length(x.win.(CHR{j})(1,:));
    l = 1;
    
    for k = 1:length(A)
        if k > 1
            if abs(A(k)) < abs(A(k-l))
                if A(k) == 0
                    l = l + 1;
                    continue
                else
                    j = j + 1;
                    len = length(x.win.(CHR{j})(1,:));
                    l = 1;
                end
            elseif l > 1
                j = j + 1;
                len = length(x.win.(CHR{j})(1,:));
                l = 1;
            end
        end

        if A(k) > 0
            strand = '+';
        else
            strand = '-';
        end
        
        if A(k) == 0
            A(k) = 1;
            read_length = read_length - 1;
        end
        if strand == '+';
            win = abs(A(k)):abs(A(k)) + read_length - 1;
        else
            win = abs(A(k)) - read_length + 1:abs(A(k));
        end
        if win(1) < 1 || round(win(end)/window) > x.win.(CHR{j})(1,end)
            continue
        end
        
        x.bp.(CHR{j})(k) = win(1);
        ls = round(win(1)/window);
        
        if ls < 1
            ls = 1;
        end
        if ls > len
            x.win.(CHR{j})(2,end) = x.win.(CHR{j})(2,end) + 1;
            continue
        end
        if win(1) < x.win.(CHR{j})(1,ls)
            while win(1) < x.win.(CHR{j})(1,ls)
                ls = ls + 1;
            end
            if win(end) > x.win.(CHR{j})(1,ls)
                x.win.(CHR{j})(2,ls) = x.win.(CHR{j})(2,ls) + 1;
            else
                a = x.win.(CHR{j})(1,ls+1) - win(1);
                b = a/read_length;
                x.win.(CHR{j})(2,ls) = x.win.(CHR{j})(2,ls) + b;
                x.win.(CHR{j})(2,ls+1) = x.win.(CHR{j})(2,ls+1) + 1 - b;
            end
        else
            if win(1) >= x.win.(CHR{j})(1,end)
                x.win.(CHR{j})(2,end) = x.win.(CHR{j})(2,end) + 1;
                continue
            end
            while win(1) >= x.win.(CHR{j})(1,ls)
                ls = ls + 1;
            end
            if win(end) < x.win.(CHR{j})(1,ls)
                x.win.(CHR{j})(2,ls) = x.win.(CHR{j})(2,ls) + 1;
            else
                a = win(end)-x.win.(CHR{j})(1,ls)+1;
                b = a/read_length;
                x.win.(CHR{j})(2,ls-1) = x.win.(CHR{j})(2,ls-1) + 1 - b;
                x.win.(CHR{j})(2,ls) = x.win.(CHR{j})(2,ls) + b;
            end
        end
        
        if mod(k,100000) == 0
            fprintf('Read %d tags\n',k)
        end
    end
    for j = 1:length(CHR)
        x.bp.(CHR{j}) = x.bp.(CHR{j})(x.bp.(CHR{j}) ~= 0);
    end
    
    
    assignin('base',headers{i},x)
end

TagCount
save(sprintf('~/ChIP-Data/%s/%s_FullGenome_%s_PostSPP',person{1},perID{1},date),headers{:},'headers','len', ...
    'read_length','window','ReadSum','test','perID','-v7.3')

for k = 1:length(headers)
    save(sprintf('~/ChIP-Data/%s/%s_%s_SumReads_PostSPP.mm9.mat',person{1},perID{1},headers{k}),sprintf('%s',headers{k}),'-v7.3')
end

clear fid ls a b x CHR chr i j tline p strand st seq phred score person A k
