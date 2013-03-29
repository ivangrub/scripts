function SumReads_Function_MacPro(aa,bb,cc,dd)

%Calculate the read coverage across the whole genome for single bp and
%condensed windows


%% Add options and identify people
person = cell(aa);
headers = cell(bb);
test = strrep(headers,'_',' ');
if ischar(cc)
    window = str2double(cc);
else
    window = cc;
end

perID = dd;


clear aa bb cc dd
data = importdata('mm9.chr.length');
len = data.data;
chr = data.textdata;

%% Initiate matrices
for i = 1:length(headers)
    for j = 1:length(chr)
        bp.(chr{j}) = zeros(1,round(len(j)./window));
        left = 1:window:len(j)-window;
        right = zeros(1,round(len(j)./window));
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
    assignin('caller',headers{i},Y);
    clear bp win Y
end

%% Scan through the data
for i = 1:length(headers)
    fprintf(sprintf('Working on %s\n',headers{i}))
    fid = fopen(sprintf('/Volumes/Ivan/Genomic/ChIP-Seq/%s/%s_%s.mm9.bw',person{1},person{2},headers{i}),'r');
    x = evalin('caller',headers{i});
    j = 0;
    while 1
        j = j + 1;
        tline = fgetl(fid);
        if ~ischar(tline)
            break
        end
        [~,strand,CHR,st,seq,~,~] = strread(tline,'%s%c%s%d%s%s%d','delimiter','\t');
        read_length = length(seq{1});
        if st == 0
            st = 1;
            read_length = read_length - 1;
        end
        if strand == '+';
            win = st:st + read_length - 1;
        else
            win = st - read_length + 1:st;
        end
        if min(win) < 1 || max(round(win./window)) > x.win.(CHR{1})(1,end)
            continue
        end
        x.bp.(CHR{1})(j) = win(1);
        ls = round(win(1)./window);
        if ls < 1
            ls = 1;
        end
        if win(1) < x.win.(CHR{1})(1,ls)
            while win(1) < x.win.(CHR{1})(1,ls)
                ls = ls + 1;
            end
            if win(end) > x.win.(CHR{1})(1,ls)
                x.win.(CHR{1})(2,ls) = x.win.(CHR{1})(2,ls) + 1;
            else
                a = x.win.(CHR{1})(1,ls+1) - win(1);
                b = a/read_length;
                x.win.(CHR{1})(2,ls) = x.win.(CHR{1})(2,ls) + b;
                x.win.(CHR{1})(2,ls+1) = x.win.(CHR{1})(2,ls+1) + 1 - b;
            end
        else
            if win(1) >= x.win.(CHR{1})(1,end)
                x.win.(CHR{1})(2,end) = x.win.(CHR{1})(2,end) + 1;
                continue
            end
            while win(1) >= x.win.(CHR{1})(1,ls)
                ls = ls + 1;
            end
            if win(end) < x.win.(CHR{1})(1,ls)
                x.win.(CHR{1})(2,ls) = x.win.(CHR{1})(2,ls) + 1;
            else
                a = win(end)-x.win.(CHR{1})(1,ls)+1;
                b = a./read_length;
                x.win.(CHR{1})(2,ls-1) = x.win.(CHR{1})(2,ls-1) + 1 - b;
                x.win.(CHR{1})(2,ls) = x.win.(CHR{1})(2,ls) + b;
            end
        end
        
        if mod(j,100000) == 0
            fprintf('Read %d tags\n',j)
        end
    end
    for j = 1:35
        x.bp.(chr{j}) = x.bp.(chr{j})(x.bp.(chr{j}) ~= 0);
    end
    
    
    eval([genvarname(headers{i}) '= x;']);
    fclose(fid);
end
clear x

ReadSum = zeros(1,length(headers));
for i = 1:length(headers)
    x = eval(headers{i});
    chr = fieldnames(x.win);
    for j = 1:length(chr)
        ReadSum(i) = ReadSum(i) + sum(x.win.(chr{j})(2,:));
    end
end

save(sprintf('/Volumes/Ivan/Genomic/ChIP-Seq/%s_FullGenome_%s',perID{1},date),headers{:},'headers','len', ...
    'read_length','window','ReadSum','test','perID','-v7.3')

for k = 1:length(headers)
    save(sprintf('/Volumes/Ivan/Genomic/ChIP-Seq/%s_%s_SumReads.mm9.mat',perID{1},headers{k}),sprintf('%s',headers{k}))
end

clear fid ls a b x CHR chr i j tline p strand st seq phred score person
