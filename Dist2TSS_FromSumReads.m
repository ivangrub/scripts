%Find the Distance to nearest TSS


x = eval(headers{1});

%Find the minimum distance to the closest TSS for each of the enriched
%intergenic regions with respect to each IP

chr = fieldnames(x.win);

fprintf('Initialize Dist structure\n')
conden_len = zeros(35,1);
Dist = NaN(length(x.win.chr1),10);
Chr = ones(1,length(x.win.chr1(1,:)));
Win = x.win.chr1(1,:)+window/2;
conden_len(1) = length(x.win.chr1);
for j = 2:length(chr)
    Dist = [Dist;NaN(length(x.win.(chr{j})),10)];
    Win = [Win x.win.(chr{j})(1,:)+read_length/2];
    Chr = [Chr j*ones(1,length(x.win.(chr{j})))];
    conden_len(j) = length(x.win.(chr{j}));
end

l = 0;
k = 1;

for j = 1:length(chr)
    fprintf('On %s\n',chr{j});
    chr_TSS = strcmp(chr{j},knownGene(:,4));
    o = find(chr_TSS);
    O = length(o);
    TSS = zeros(O,1);
    TES = zeros(O,1);
    pc = find(cell2mat(knownGene(o,5)) == 1);
    nc = find(cell2mat(knownGene(o,5)) == -1);
    TSS(pc) = cell2mat(knownGene(o(pc),6));
    TES(pc) = cell2mat(knownGene(o(pc),7));
    TSS(nc) = cell2mat(knownGene(o(nc),7));
    TES(nc) = cell2mat(knownGene(o(nc),6));
    for i = 1:sum(conden_len(j))
        
        %Distance to TSS
        Dist2Win = Win(1,k) - TSS;
        left = (sign(Dist2Win) == -1);                      %Identify those upstream
        right = (sign(Dist2Win) == 1);                      %Identify those downstream
        [left_dist2tss,~] = max(Dist2Win(left));            %Find closest upstream
        [right_dist2tss,~] = min(Dist2Win(right));          %Find closest downstream
        if isempty(left_dist2tss) == 0 && isempty(right_dist2tss) == 0
            geneup = (left_dist2tss == Win(1,k) - TSS);
            b = find(geneup,1,'first');
            geneleft = l + b;
            genedown = (right_dist2tss == Win(1,k) - TSS);
            c = find(genedown,1,'first');
            generight = l + c;
            Dist(k,1:5) = [Win(1,k) geneleft left_dist2tss generight right_dist2tss];
        elseif isempty(left_dist2tss) == 0 && isempty(right_dist2tss) == 1
            geneup = (left_dist2tss == Win(1,k) - TSS);
            b = find(geneup,1,'first');
            geneleft = l + b;
            Dist(k,1:3) = [Win(1,k) geneleft left_dist2tss];
        else
            genedown = (right_dist2tss == Win(1,k) - TSS);
            c = find(genedown,1,'first');
            generight = l + c;
            Dist(k,[1 4 5]) = [Win(1,k) generight right_dist2tss];
        end
        
        
        %Distance to TES
        Dist2WinTES = Win(1,k) - TES;
        leftTES = (sign(Dist2WinTES) == -1);                %Identify those upstream
        rightTES = (sign(Dist2WinTES) == 1);                %Identify those downstream
        [left_dist2tes,~] = max(Dist2WinTES(leftTES));      %Find closest upstream
        [right_dist2tes,~] = min(Dist2WinTES(rightTES));    %Find closest downstream
        if isempty(left_dist2tes) == 0 && isempty(right_dist2tes) == 0
            geneuptes = (left_dist2tes == Win(1,k) - TES);
            btes = find(geneuptes,1,'first');
            genelefttes = l + btes;
            genedowntes = (right_dist2tes == Win(1,k) - TES);
            ctes = find(genedowntes,1,'first');
            generighttes = l + ctes;
            Dist(k,6:10) = [Win(1,k) genelefttes left_dist2tes generighttes right_dist2tes];
        elseif isempty(left_dist2tes) == 0 && isempty(right_dist2tes) == 1
            geneuptes = (left_dist2tes == Win(1,k) - TES);
            btes = find(geneuptes,1,'first');
            genelefttes = l + btes;
            Dist(k,6:8) = [Win(1,k) genelefttes left_dist2tes];
        else
            genedowntes = (right_dist2tes == Win(1,k) - TES);
            ctes = find(genedowntes,1,'first');
            generighttes = l + ctes;
            Dist(k,[6 9 10]) = [Win(1,k) generighttes right_dist2tes];
        end
        k = k + 1;
    end
    l = l + O;
end


clear i k m j o gene TSS right_dist2tss left_dist2tss c b chr chr_TSS ...
    Dist2Win Win left right find_left find_right l pc nc TES Dist2WinTES ...
    leftTES rightTES find_leftTES find_rightTES left_dist2tes btes ctes ...
    right_dist2tes O geneleft generight geneleftTES generightTES geneup ...
    geneuptes genedown genedowntes genelefttes generighttes