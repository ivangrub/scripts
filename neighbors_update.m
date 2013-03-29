%Identify neighboring windows that are enriched at a statistically
%significant level

function New = neighbors_update(len,window_thresh,x,knownGene,PI)

Conden.Distance = NaN(4*len,10);
Conden.Window = NaN(4*len,2);
Conden.Reads = NaN(4*len,1);
Conden.Enrichment = NaN(4*len,1);

a = 1;
columns = [1 1 6 6];
TSSMatrix = [3 5]; 
TESMatrix = [8 10];
for i = 1:length(knownGene)
    [X,Y] = find(x.Distance(:,[2 4 7 9]) == i);
    strand = sign(knownGene{i,5});
    
    if strand == 1
        start = 6;
        endit = 7;
    else 
        start = 7;
        endit = 6;
    end
    
    r = 1;
    for j = 1:length(X)
        dist = x.Distance(X,columns(Y)) - x.Distance(X(r),columns(Y(r)));
        if length(dist) == 1
            if Y(r) == 1 || Y(r) == 2
                Conden.Distance(a,[1 TSSMatrix(Y(r))-1 TSSMatrix(Y(r))]) = x.Distance(X(r),[1 TSSMatrix(Y(r))-1 TSSMatrix(Y(r))]);
            else
                Conden.Distance(a,[6 TESMatrix(Y(r))-1 TESMatrix(Y(r))]) = x.Distance(X(r),[6 TESMatrix(Y(r))-1 TESMatrix(Y(r))]);
            end
            Conden.Reads(a,1) = x.Reads(X(r),1);
            Conden.Window(a,1:2) = x.Window(X(r),1:2);
            Conden.Enrichment(a,1) = x.Enrichment(X(r),1);
            a = a + 1;
            break
        else
            [a1,b1] = find((abs(dist) <= window_thresh & abs(dist) > 0),1,'first');
            [a2,b2] = find((abs(dist) <= window_thresh & abs(dist) > 0),1,'last');
            if isempty(a1) == 0 && isempty(a2) == 0
                if a1 == a2 && b1 == b2
                    a1 = [];
                    if Y(r) == 1 || Y(r) == 2
                        Conden.Distance(a,[1 TSSMatrix(Y(a2))-1 TSSMatrix(Y(b2))]) = x.Distance(X(a2),[1 TSSMatrix(Y(b2))-1 TSSMatrix(Y(a2))]);
                    else
                        Conden.Distance(a,[6 TESMatrix(Y(a2))-1 TESMatrix(Y(b2))]) = x.Distance(X(a2),[6 TESMatrix(Y(b2))-1 TESMatrix(Y(a2))]);
                    end
                    Conden.Distance(a,:) = x.Distance(dist(a2),DistMatrix(Y(a2)));
                    Conden.Reads(a,1) = x.Reads(dist(a2),:);
                    Conden.Window(a,1:2) = x.Window(dist(a2),1:2);
                    Conden.Enrichment(a,1) = x.Enrichment(dist(a2),1);
                    a = a + 1;
                end
                if isempty(a1) == 0
                    Conden.Window(a,1:2) = [x.Window(dist(a1),1) x.Window(dist(a2),2)];
                    if Y(r) == 1 || Y(r) == 2
                        Conden.Distance(a,[1 TSSMatrix(Y(a2))-1 TSSMatrix(Y(a2))]) = [mean(Conden.Window(a,1:2)) i (mean(Conden.Window(a,1:2)) - knownGene{i,start})];
                    else
                        Conden.Distance(a,[6 TESMatrix(Y(a2))-1 TESMatrix(Y(a2))]) = [mean(Conden.Window(a,1:2)) i (mean(Conden.Window(a,1:2)) - knownGene{i,endit})];
                    end
                    if j == 1 || a == 1
                        Conden.Reads(a,1) = sum(x.Reads(dist(a1:a2),1));
                        Conden.Enrichment(a,1) = Conden.Reads(a,1)/sum(PI(dist(a1:a2),1));
                        a = a + 1;
                    elseif j > 1 && (Conden.Window(a,1) >= Conden.Window(a-1,2))
                        Conden.Reads(a,1) = sum(x.Reads(dist(a1:a2),1));
                        Conden.Enrichment(a,1) = Conden.Reads(a,1)/sum(PI(dist(a1:a2),1));
                        a = a + 1;
                    else
                        Conden.Window(a,:) = NaN;
                        Conden.Distance(a,:) = NaN;
                    end
                end
            end
        end
        if a2 == length(X)
            break
        end
        r = a2 + 1;
    end
end

[A,B] = find(isnan(Conden.Distance) == 0);
Conden.Distance = Conden.Distance(A,B);
Conden.Reads = Conden.Reads(A);
Conden.Window = Conden.Window(A,:);
Conden.Enrichment = Conden.Enrichment(A);

New = Conden;
end

