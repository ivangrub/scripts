function Conden = neighbors_fromwindows_v2(len,window_thresh,win,x,knownGene,PI)

Conden = NaN(4*len,8);

edge = floor(window_thresh/win);
Mid = mean(x.Window,2);

a = 1;
b = 1;
for r = 1:len-1  
    
    if (r + edge) < len
        dist = x.Window(r:r+edge,1) - x.Window(r,2);
        l = r:r+edge;
    else 
        dist = x.Window(r:end,1) - x.Window(r,2);
        l = r:len;
    end
    a1 = find(abs(dist) <= window_thresh-1,1,'first');
    a2 = find(abs(dist) <= window_thresh-1,1,'last');
    a1 = l(a1);
    a2 = l(a2);
    
    %For when it is a single window
    if a1 == a2
        Conden(a,1:2) = x.Window(a1,1:2);
        if a > 1
            if Conden(a,2) == Conden(a-1,2)
                b = b + 1;
                continue
            end
        end
        if b > 1 && (Conden(a,1) - Conden(a-1,2) <= 0) && (abs(Conden(a,1) - Conden(a-1,2)) <= window_thresh)
            b = b - 1;
            A = unique(x.Distance(b:a2,[2 4 7 9]));
            l = isnan(A) == 0;
            A = A(l);
            B = unique(x.Distance(b,[2 4 7 9]));
            l = isnan(B) == 0;
            l = sum(l);
            a = a - l - 1;
            if isempty(A) == 0
                for k = 1:length(A)
                    Conden(a,1:2) = [x.Window(b,1) x.Window(a1,2)];
                    Conden(a,7) = sum(x.Reads(1,b:a2));
                    Conden(a,8) = Conden(a,7)./sum(PI(1,b:a2));
                    NewMid = mean(Conden(a,1:2),2);
                    if knownGene{A(k),5} == 1
                        start = 6;
                        endit = 7;
                    else
                        start = 7;
                        endit = 6;
                    end
                    TSS = NewMid - knownGene{A(k),start};
                    TES = NewMid - knownGene{A(k),endit};
                    Conden(a,3:6) = [NewMid A(k) TSS TES];
                    a = a + 1;
                end
            end
            b = b + 2;
            continue
        else
            A = unique(x.Distance(a1,[2 4 7 9]));
            l = isnan(A) == 0;
            A = A(l);
            if isempty(A) == 0
                for k = 1:length(A)
                    Conden(a,1:2) = x.Window(a1,1:2);
                    Conden(a,7) = x.Reads(a1,1);
                    Conden(a,8) = x.Enrichment(a1,1);
                    if knownGene{A(k),5} == 1
                        start = 6;
                        endit = 7;
                    else
                        start = 7;
                        endit = 6;
                    end
                    TSS = Mid(a1) - knownGene{A(k),start};
                    TES = Mid(a1) - knownGene{A(k),endit};
                    Conden(a,3:6) = [Mid(a1) A(k) TSS TES];
                    a = a + 1;
                end
            end
            b = b + 1;
            continue
        end
    end
    
    %Condense the windows
    Conden(a,1:2) = [x.Window(a1,1) x.Window(a2,2)];
    if a > 1
        if Conden(a,2) == Conden(a-1,2)
            b = b + 1;
            continue
        end 
    end
    if b > 1 && (Conden(a,1) - Conden(a-1,2) <= 0) && (abs(Conden(a,1) - Conden(a-1,2)) <= window_thresh)
        b = b - 1;
        A = unique(x.Distance(b:a2,[2 4 7 9]));
        l = isnan(A) == 0;
        A = A(l);
        B = unique(x.Distance(b,[2 4 7 9]));
        l = isnan(B) == 0;
        l = sum(l);
        a = a - l - 1;
        if isempty(A) == 0
            for k = 1:length(A)
                Conden(a,1:2) = [x.Window(b,1) x.Window(a2,2)];
                Conden(a,7) = sum(x.Reads(b:a2,1));
                Conden(a,8) = Conden(a,7)./sum(PI(b:a2,1));
                NewMid = mean(Conden(a,1:2),2);
                if knownGene{A(k),5} == 1
                    start = 6;
                    endit = 7;
                else
                    start = 7;
                    endit = 6;
                end
                TSS = NewMid - knownGene{A(k),start};
                TES = NewMid - knownGene{A(k),endit};
                Conden(a,3:6) = [NewMid A(k) TSS TES];
                a = a + 1;
            end
        end
        b = b + 2;
        continue
    else
        Conden(a,7) = sum(x.Reads(a1:a2,1));
        Conden(a,8) = Conden(a,7)./sum(PI(a1:a2,1));
        NewMid = mean(Conden(a,1:2),2);
        A = unique(x.Distance(a1:a2,[2 4 7 9]));
        l = isnan(A) == 0;
        A = A(l);
        if isempty(A) == 0
            for k = 1:length(A)
                Conden(a,1:2) = [x.Window(a1,1) x.Window(a2,2)];
                Conden(a,7) = sum(x.Reads(a1:a2,1));
                Conden(a,8) = Conden(a,7)./sum(PI(a1:a2,1));
                if knownGene{A(k),5} == 1
                    start = 6;
                    endit = 7;
                else
                    start = 7;
                    endit = 6;
                end
                TSS = NewMid - knownGene{A(k),start};
                TES = NewMid - knownGene{A(k),endit};
                Conden(a,3:6) = [NewMid A(k) TSS TES];
                a = a + 1;
            end
        end
        b = b +1;
    end
end

A = (isnan(Conden(:,3)) == 0);
Conden = Conden(A,:);
    
end