clear

%Provide details about the data that needs to be inputted
%Cell Type and then all of the subsequent factors that were tested

options = {'ESC','TBP','PI'};

%Compile known gene transcription start sites (TSS)
x = compile_known;

%Add the tRNA annotations to the end of the known matrix

j = length(x);
A = importdata('tRNAannot.csv',',');
for i = 1:length(A)
    yy = strread(A{i},'%s','delimiter',',');
    x{j+i,2} = yy{2};
    x{j+i,4} = strrep(yy{1},'chr','');
    x{j+i,6} = str2double(yy{3});
    x{j+i,7} = str2double(yy{4});
    x(j+i,[1 3 5 8 9 10 11]) = {NaN NaN NaN NaN NaN NaN NaN};
end

%Add the CpG data to the knownGene data
A = importdata('cpgIslandExt.csv',',');
textdata = A.textdata;
data = A.data;
[~,indice] = sort(data,1,'descend');
data = data(indice,:);
textdata = textdata(indice,:);
CpGfrom = floor(str2double(textdata(:,3))/25);
CpGto = floor(str2double(textdata(:,4))/25);
CpGchr = strrep(textdata(:,2),'chr','');
CpG = zeros([length(x),2]);
for i = 1:length(x)
    tss = floor(cell2mat(x(i,6))/25);
    chr = cell2mat(x(i,4)); 
    I = (strcmp(CpGchr,chr) & (tss - 20) <= CpGto & (tss + 20) >= CpGfrom);
    if any(I) 
        j = find(I,1,'first');
        CpG(i,1:2) = [data(j,2) data(j,6)];
    end
end

%Iterate through the different factors (TBP, TAF1, PolII_Ser5, etc.)
n = 1;
for i = 2:length(options)
    clear data
    type = cell2mat(options(1));    
    f = cell2mat(options(i));
    data = open(sprintf('Francisco/%s_%s.mm9.mat',type,f));
    
    %Initialize while loop which will look at each TSS and determine the 
    %enrichment
    j = 1;
    next = 0;
    while next == 0
        tss = cell2mat(x(j,6));
        chr = cell2mat(x(j,4));
        CellType(j,n:n+1) = coverage(chr,tss,data);
        if j == length(x)
            next = 1;
        end
        j = j+1;
    end
    n = n + 2;
end
knownGene = x;
%exportdata(knownGene,CpG,CellType,options{1})
