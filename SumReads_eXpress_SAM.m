%Calculate the read coverage across the whole genome for single bp and
%condensed windows


%% Add options and identify people
person = {'Claudia' 'CC'};
headers = {'Rad23b' 'PI'};
window = 250;

data = importdata('mm9.chr.length');
len = data.data;
chr = data.textdata;

%% Initiate matrices
for i = 1:length(headers)
    for j = 1:length(chr)
        bp.(chr{j}) = zeros(2,round(len(j)/window));
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

%% Scan through the data
for i = 1:length(headers)
    fprintf(sprintf('Working on %s\n',headers{i}))
    fid = fopen(sprintf('%s/%s_%s.mm9.bw',person{1},person{2},headers{i}),'r');
    x = eval(headers{i});
    j = 0;
    while 1
        j  = j + 1;
        tline = fgetl(fid);
        if ~ischar(tline)
            break
        end
        [tag,~,bin,st,~,~,~,~,~,~,~,seq,~,~,frac_str] = strread(tline,'%s%s%s%d%s%s%s%s%s%s%s%s%s%s%s','delimiter','\t');
        if ~strcmp(tag{1}[1],'H')
        	continue
        end
        ch = textscan(bin,'delimiter','!');
        CHR = ch{2};
        read_length = length(seq);
        frac = str2num(strrep('XP:f:',frac_str)) ;
        if st == 0
            st = 1;
            read_length = read_length - 1;
        end
        win = st:st + read_length - 1;
        if min(win) < 1 || max(round(win/window)) > x.win.(CHR)(1,end)
             continue
        end
        x.bp.(CHR)(1:2,j) = [win(1);frac];
        ls = round(win(1)/window);
        if ls < 1
            ls = 1;
        end
        if win(1) < x.win.(CHR)(1,ls)
            while win(1) < x.win.(CHR)(1,ls) 
                ls = ls + 1;
            end
            if win(end) > x.win.(CHR)(1,ls)
                x.win.(CHR)(2,ls) = x.win.(CHR)(2,ls) + frac;
            else
                a = x.win.(CHR)(1,ls+1) - win(1);
                b = a/read_length;
                x.win.(CHR)(2,ls) = x.win.(CHR)(2,ls) + b*frac;
                x.win.(CHR)(2,ls+1) = x.win.(CHR)(2,ls+1) + frac - b*frac;
            end
        else
            if win(1) >= x.win.(CHR)(1,end)
                x.win.(CHR)(2,end) = x.win.(CHR)(2,end) + frac;
                continue
            end
            while win(1) >= x.win.(CHR)(1,ls)
                ls = ls + 1;
            end
            if win(end) < x.win.(CHR)(1,ls)
                x.win.(CHR)(2,ls) = x.win.(CHR)(2,ls) + frac;
            else
                a = win(end)-x.win.(CHR)(1,ls)+1;
                b = a/read_length;
                x.win.(CHR)(2,ls-1) = x.win.(CHR)(2,ls-1) + frac - b*frac;
                x.win.(CHR)(2,ls) = x.win.(CHR)(2,ls) + b*frac;
            end
        end 
        
        if mod(j,100000) == 0
            fprintf('Read %d tags\n',j)
        end
    end
    for j = 1:35
        x.bp.(chr{j}) = x.bp.(chr{j})(x.bp.(chr{j}) ~= 0);
    end
    assignin('base',headers{i},x)
    fclose(fid);
end
TagCount
save(sprintf('%s_FullGenome_%s',person{2},date),headers{:},'headers','len', ...
    'read_length','window','ReadSum','-v7.3')
for i = 1:length(headers)
	save(sprintf('%s_SumReads.mm9.mat',headers{i}),'-v7.3')
end
clear fid ls a b x CHR chr i j tline p strand st seq phred score person
