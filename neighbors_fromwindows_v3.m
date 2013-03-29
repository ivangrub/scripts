function [Conden Chr_conden TSS_conden] = neighbors_fromwindows_v3(len,window_thresh,x,knownGene,PI,chr)

Conden = NaN(4*len,10);
Chr_conden = NaN(1,4*len);
TSS_conden = NaN(2,4*len);

a = 1;
for r = 1:length(chr)
    fprintf('Condensing on %s\n',chr{r})
    subwin = find(x.Chr(:) == r);
   
    if isempty(subwin)
        continue
    end
    
    for f = 1:length(subwin)
        
        dist = x.Window(2,subwin) - x.Window(1,subwin(f));
        
        a1 = find(abs(dist) <= window_thresh,1,'first');
        a2 = find(abs(dist) <= window_thresh,1,'last');
        
        if isempty(a1) && isempty(a2)
            continue
        end
        
        if a1 == a2
            Conden(a,1:2) = x.Window(1:2,subwin(a1));
            
            %Addresses the issue if the last window of a sequence has already
            %been accounted for
            if a > 1 && (Conden(a,2) == Conden(a-1,2))
                continue
            end
            
            %Addresses the issue if there is an overlap in windows
            if a > 1 && (Conden(a,1) <= Conden(a-1,2)) && (Conden(a,2) >= Conden(a-1,2))
                b = find(x.Window(1,subwin) == Conden(a-1,1));
                A = unique(x.Distance(subwin(b):subwin(a1),[2 4 7 9]));
                l = isnan(A) == 0;
                A = A(l);
                l = sum(l);
                a = a - l;
                if isempty(A) == 0
                    for k = 1:length(A)
                        Conden(a,1:2) = [x.Window(1,subwin(b)) x.Window(2,subwin(a1))];
                        [Conden(a,7),summit] = max(x.Reads(1,subwin(b):subwin(a1)));
                        Conden(a,8) = sum(x.Reads(1,subwin(b):subwin(a1)))./sum(PI(1,subwin(b):subwin(a1)));
                        space = subwin(b):subwin(a2);
                        Conden(a,9:10) = [x.Window(1,space(summit)) x.Window(2,space(summit))];
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
                        TSS_conden(:,a) = [cell2mat(knownGene(A(k),start));cell2mat(knownGene(A(k),endit))];
                        Chr_conden(1,a) = r;
                        a = a + 1;
                    end
                end
                continue
            else
                A = unique(x.Distance(subwin(a1),[2 4 7 9]));
                l = isnan(A) == 0;
                A = A(l);
                if isempty(A) == 0
                    for k = 1:length(A)
                        Conden(a,1:2) = x.Window(1:2,subwin(a1));
                        Conden(a,7) = x.Reads(1,subwin(a1));
                        Conden(a,8) = x.Enrichment(1,subwin(a1));
                        Conden(a,9:10) = Conden(a,1:2);
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
                        TSS_conden(:,a) = [cell2mat(knownGene(A(k),start));cell2mat(knownGene(A(k),endit))];
                        Chr_conden(1,a) = r;
                        a = a + 1;
                    end
                end
                continue
            end
        end
        
        %If there are windows to condense
        Conden(a,1:2) = [x.Window(1,subwin(a1)) x.Window(2,subwin(a2))];
        
        %Addresses the issue if the last window of a sequence has already
        %been accounted for
        if a > 1 && (Conden(a,2) == Conden(a-1,2))
            continue
        end
        
        %Addresses the issue if there is an overlap in windows
        if a > 1 && (Conden(a,1) <= Conden(a-1,2)) && (Conden(a,2) >= Conden(a-1,2))
            b = find(x.Window(1,subwin) == Conden(a-1,1));
            A = unique(x.Distance(subwin(b):subwin(a2-1),[2 4 7 9]));
            l = isnan(A) == 0;
            A = A(l);
            l = sum(l);
            a = a - l;
            if isempty(A) == 0
                for k = 1:length(A)
                    Conden(a,1:2) = [x.Window(1,subwin(b)) x.Window(2,subwin(a2))];
                    [Conden(a,7),summit] = max(x.Reads(1,subwin(b):subwin(a2)));
                    Conden(a,8) = sum(x.Reads(1,subwin(b):subwin(a2)))./sum(PI(1,subwin(b):subwin(a2)));
                    space = subwin(b):subwin(a2);
                    Conden(a,9:10) = [x.Window(1,space(summit)) x.Window(2,space(summit))];
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
                    TSS_conden(:,a) = [cell2mat(knownGene(A(k),start));cell2mat(knownGene(A(k),endit))];
                    Chr_conden(1,a) = r;
                    a = a + 1;
                end
            end
            continue
        else
            A = unique(x.Distance(subwin(a1):subwin(a2),[2 4 7 9]));
            l = isnan(A) == 0;
            A = A(l);
            if isempty(A) == 0
                for k = 1:length(A)
                    Conden(a,1:2) = [x.Window(1,subwin(a1)) x.Window(2,subwin(a2))];
                    [Conden(a,7),summit] = max(x.Reads(1,subwin(a1):subwin(a2)));
                    Conden(a,8) = sum(x.Reads(1,subwin(a1):subwin(a2)))./sum(PI(1,subwin(a1):subwin(a2)));
                    space = subwin(a1):subwin(a2);
                    Conden(a,9:10) = [x.Window(1,space(summit)) x.Window(2,space(summit))];
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
                    TSS_conden(:,a) = [cell2mat(knownGene(A(k),start));cell2mat(knownGene(A(k),endit))];
                    Chr_conden(1,a) = r;
                    a = a + 1;
                end
            end
        end
    end
end

A = (isnan(Conden(:,3)) == 0);
Conden = Conden(A,:);
infini = Conden(:,8) == Inf;
Conden(infini,8) = Conden(infini,7);

A = isnan(Chr_conden(1,:)) == 0;
Chr_conden = Chr_conden(1,A);
end