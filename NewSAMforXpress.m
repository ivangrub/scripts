IP = 'MEF_TBP';
person = {'Teppei' 'TY'};

bin = 1000;
overlap = 60;

data = importdata('mm9.chr.length');
len = data.data;
chr = data.textdata;


%% Initiate matrices

for j = 1:length(chr)
    left = 1:bin-overlap:len(j)-bin-overlap;
    right = bin:bin-overlap:len(j);
    if length(left) ~= length(right)
        if length(left) > length(right)
            left = left(1:length(right));
        else
            left(length(left)+1) = left(end) + bin + overlap;
        end
    end
    value = zeros(1,length(left));
    win.(chr{j}) = [left;right;value];
end

%% Process SAM and FASTA files
fid = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s.mm9.sam',person{1},person{2},IP),'r');
fid2 = fopen(sprintf('eXpress_%s.mm9.sam',IP),'w');
fid3 = fopen(sprintf('eXpress_%s.mm9.fasta',IP),'w');

j = 0;
z = 1;
while 1
    j = j + 1;
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    if j <= 24
       fprintf(fid2,'%s\n',tline); 
       continue
    end
    
    A = strread(tline,'%s','delimiter','\t');
    
    if mod(j,100000) == 0
        fprintf(sprintf('%d Reads Processed\n',j))
    end
    
    CHR = A{3,1};
    if strcmp(CHR,'*')
        continue
    end
    st = str2double(A{4,1});
    
    [~,pos] = min(abs(win.(CHR)(1,:) - st));
    
    if pos > 1
        if (st - win.(CHR)(1,pos)) <= 0 
            a = pos - 1;
        elseif (st > win.(CHR)(1,pos)) && (st < win.(CHR)(2,pos-1))
            a = [pos - 1,pos];
        else
            a = pos;
        end
    else
        a = pos;
    end
    
    for k = 1:length(a)
        %Create ID name
        abs_bin = sprintf('bin%d',z);
        position = strcat(abs_bin,sprintf('!%s!%d!%d',CHR,win.(CHR)(1,a(k)),win.(CHR)(2,a(k))));
        offset = st - win.(CHR)(1,a(k));
        for Y = 1:2
            fprintf(fid2,'%s\t',A{Y});
        end
        fprintf(fid2,'%s\t%d\t',position,offset);
        for Y = 5:length(A)-1
            fprintf(fid2,'%s\t',A{Y});
        end  
        fprintf(fid2,'%s\n',A{Y+1});
        z = z + 1;
        
        %Determine that the bin has reads
        win.(CHR)(3,a(k)) = 1;
        
    end
    
end

fprintf('Printing FASTA File')
for i = 1:length(chr)
    x = fastaread('MM9CHR/%s.fa',chr{i});
    L = win.(chr{i})(:,win.(chr{i})(3,:)==1);
    for j = 1:sum(L)
        index = round(linspace(L(1,j),L(2,j),20));
        fprintf(fid3,sprintf('>%s\n%s\n',position,upper(x.Sequence(index(1):index(2)))));
        for C = 2:length(index)-1
            fprintf(fid3,sprintf('%s\n',upper(x.Sequence(index(C)+1:index(C+1)))));
        end
        fprintf(fid3,'\n');
    end
end
fclose(fid);fclose(fid2);fclose(fid3);




