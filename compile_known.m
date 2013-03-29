function known = compile_known

A = importdata('knownGene.csv',',');
data = A.data;
textdata = A.textdata;
B = importdata('kgXref.csv',',');

[m,n] = size(data);
known(:,1) = textdata(:,1);
accnum = cellstr(known(:,1));
for i = 1:m
    %Input the name of the gene
    yy = strread(B{i},'%s','delimiter',',');
    I = strcmp(accnum,cellstr(yy(1)));
    j = find(I,1,'first');
    known(j,2) = cellstr(yy(5));
    known(j,3) = strrep(cellstr(yy(8)),'"','');
    
    %Convert Chr# to a string reflecting the chromosome number
    known(i,4) = strrep(textdata(i,2),'chr','');
    
    %Convert the plus/minus to 1 or -1 respectively
    if strcmp(textdata(i,3),'+') == 1
        known(i,5) = num2cell(1);
    else
        known(i,5) = num2cell(-1);
    end
end

%Inputs the TSS, TSE, CSS and CES numbers 
known(:,6:9) = num2cell(data);

for i = 1:m-2
    %Is the TSS unique? yes/no = 1/0
    if cell2mat(known(i,6)) == cell2mat(known(i+1,6))
        known(i,10) = num2cell(0);
    else
        known(i,10) = num2cell(1);
    end
    if i > 1
       if cell2mat(known(i,6)) == cell2mat(known(i-1,6))
           known(i,10) = num2cell(0);
       else
           known(i,10) = num2cell(1);
       end
       if cell2mat(known(i,10)) ~= cell2mat(known(i-1,10))
            if cell2mat(known(i,6)) == cell2mat(known(i-1,6))
                known(i-1,10) = num2cell(0);
            end
       end
    end
    %Is the TES unique? yes/no = 1/0
    if cell2mat(known(i,7)) == cell2mat(known(i+1,7))   
        known(i,11) = num2cell(0);
    else
        known(i,11) = num2cell(1);
    end
    if i > 1
       if cell2mat(known(i,7)) == cell2mat(known(i-1,7))
           known(i,11) = num2cell(0);
       else
           known(i,11) = num2cell(1);
       end
       if cell2mat(known(i,11)) ~= cell2mat(known(i-1,11))
            if cell2mat(known(i,7)) == cell2mat(known(i-1,7))
                known(i-1,11) = num2cell(0);
            end
       end
    end 
end
end