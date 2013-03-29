function Conden = neighbors_fromwindows(len,window_thresh,win,x,knownGene,PI)

Conden.Window = NaN(len,2);
Conden.Distance = NaN(4*len,4);
Conden.Reads = NaN(len,1);
Conden.Enrichment = NaN(len,1);

edge = floor(window_thresh/win);
Mid = mean(x.Window,2);

a = 1;
b = 1;
for r = 1:len-1  
    
    if (r + edge) < len
        dist = x.Window(r:r+edge,2) - x.Window(r,1);
        l = r:r+edge;
    else 
        dist = x.Window(r:end,2) - x.Window(r,1);
        l = r:len;
    end
    a1 = find(abs(dist) <= window_thresh,1,'first');
    a2 = find(abs(dist) <= window_thresh,1,'last');
    a1 = l(a1);
    a2 = l(a2);
    
    %For when it is a single window
    if a1 == a2
        Conden.Window(a,1:2) = x.Window(a1,1:2);
        Conden.Reads(a,1) = x.Reads(a1,1);
        Conden.Enrichment(a,1) = x.Enrichment(a1,1);
        A = unique(x.Distance(a1,[2 4 7 9]));
        l = isnan(A) == 0;
        A = A(l);
        if isempty(A) == 0
            for j = 1:length(A)
                if knownGene{A(j),5} == 1
                    start = 6;
                    endit = 7;
                else
                    start = 7;
                    endit = 6;
                end
                TSS = Mid(a1) - knownGene{A(j),start};
                TES = Mid(a1) - knownGene{A(j),endit};
                Conden.Distance(b,1:4) = [Mid(a1) A(j) TSS TES];
                b = b + 1;
            end
        end
        if a > 1 && (Conden.Window(a,1) < Conden.Window(a-1,2)) && (abs(Conden.Window(a-1,2)-Conden.Window(a,1)) <= window_thresh)
            a = a - 1;
            b = b - length(A);
            Conden.Window(a,1:2) = [Conden.Window(a,1) Conden.Window(a+1,2)];
            Conden.Reads(a,1) = sum(Conden.Reads(a:a+1,1));
            Conden.Enrichment(a,1) = Conden.Reads(a,1)./sum(PI(a:a+1,1));
            A = unique(x.Distance(a,[2 4 7 9]));
            l = isnan(A) == 0;
            A = A(l);
            if isempty(A) == 0
                for j = 1:length(A)
                    if knownGene{A(j),5} == 1
                        start = 6;
                        endit = 7;
                    else
                        start = 7;
                        endit = 6;
                    end
                    TSS = Mid(a1) - knownGene{A(j),start};
                    TES = Mid(a1) - knownGene{A(j),endit};
                    Conden.Distance(b,1:4) = [Mid(a1) A(j) TSS TES];
                    b = b + 1;
                end
            end
        end
        a = a + 1;
        continue
    end
    
    %Condense the windows
    Conden.Window(a,1:2) = [x.Window(a1,1) x.Window(a2,2)];
    Conden.Reads(a,1) = sum(x.Reads(a1:a2,1));
    Conden.Enrichment(a,1) = Conden.Reads(a,1)./sum(PI(a1:a2,1));
    NewMid = mean(Conden.Window(a,1:2),2);
    A = unique(x.Distance(a1:a2,[2 4 7 9]));
    l = isnan(A) == 0;
    A = A(l);
    if isempty(A) == 0
        for j = 1:length(A)
            if knownGene{A(j),5} == 1
                start = 6;
                endit = 7;
            else
                start = 7;
                endit = 6;
            end
            TSS = NewMid - knownGene{A(j),start};
            TES = NewMid - knownGene{A(j),endit};
            Conden.Distance(b,1:4) = [NewMid A(j) TSS TES];
            b = b + 1;
        end
    end
    if a > 1 && (Conden.Window(a,1) < Conden.Window(a-1,2)) && (abs(Conden.Window(a-1,2)-Conden.Window(a,1)) <= window_thresh)
        a = a - 1;
        b = b - length(A);
        Conden.Window(a,1:2) = [Conden.Window(a,1) Conden.Window(a+1,2)];
        Conden.Reads(a,1) = sum(Conden.Reads(a:a+1,1));
        Conden.Enrichment(a,1) = Conden.Reads(a,1)./sum(PI(a:a+1,1));
        NewMid = mean(Conden.Window(a,1:2),2);
        A = unique(x.Distance(a:a+1,[2 4 7 9]));
        l = isnan(A) == 0;
        A = A(l);
        if isempty(A) == 0
            for j = 1:length(A)
                if knownGene{A(j),5} == 1
                    start = 6;
                    endit = 7;
                else
                    start = 7;
                    endit = 6;
                end
                TSS = NewMid - knownGene{A(j),start};
                TES = NewMid - knownGene{A(j),endit};
                Conden.Distance(b,1:4) = [NewMid A(j) TSS TES];
                b = b + 1;
            end
        end
    end
    
    a = a + 1;
end

A = (isnan(Conden.Reads(:,1)) == 0);
Conden.Window = Conden.Window(A,:);
Conden.Reads = Conden.Reads(A,1);
Conden.Enrichment = Conden.Enrichment(A,1);
B = (isnan(Conden.Distance(:,1)) == 0);
Conden.Distance = Conden.Distance(B,1:4);
    
end