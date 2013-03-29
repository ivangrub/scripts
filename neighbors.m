%Identify neighboring windows that are enriched at a statistically
%significant level

function New = neighbors(len,window_thresh,x,knownGene,PI)

LeftTSS.Distance = NaN(len,3);
LeftTSS.Window = NaN(len,2);
LeftTSS.Reads = NaN(len,1);
LeftTSS.Enrichment = NaN(len,1);

RightTSS.Distance = NaN(len,3);
RightTSS.Window = NaN(len,2);
RightTSS.Reads = NaN(len,1);
RightTSS.Enrichment = NaN(len,1);

LeftTES.Distance = NaN(len,3);
LeftTES.Window = NaN(len,2);
LeftTES.Reads = NaN(len,1);
LeftTES.Enrichment = NaN(len,1);

RightTES.Distance = NaN(len,3);
RightTES.Window = NaN(len,2);
RightTES.Reads = NaN(len,1);
RightTES.Enrichment = NaN(len,1);


a = 1;
b = 1;
c = 1;
d = 1;
for i = 1:length(knownGene)
    leftTSS = find((x.Distance(:,2) == i));
    rightTSS = find((x.Distance(:,4) == i));
    leftTES = find((x.Distance(:,7) == i));
    rightTES = find((x.Distance(:,9) == i));
    strand = sign(knownGene{i,5});
    if strand == 1
        start = 6;
        endit = 7;
    else 
        start = 7;
        endit = 6;
    end
    
    r = 1;
    for j = 1:length(leftTSS)
        distleftTSS = x.Distance(leftTSS,1) - x.Distance(leftTSS(r),1);
        if length(distleftTSS) == 1
            LeftTSS.Distance(a,:) = x.Distance(leftTSS(r),1:3);
            LeftTSS.Reads(a,1) = x.Reads(leftTSS(r),1);
            LeftTSS.Window(a,1:2) = x.Window(leftTSS(r),1:2);
            LeftTSS.Enrichment(a,1) = x.Enrichment(leftTSS(r),1);
            a = a + 1;
            break
        else
            a1 = find((abs(distleftTSS) <= window_thresh),1,'first');
            a2 = find((abs(distleftTSS) <= window_thresh),1,'last');
            if a1 == a2
                a1 = [];
                LeftTSS.Distance(a,:) = x.Distance(leftTSS(a2),1:3);
                LeftTSS.Reads(a,1) = x.Reads(leftTSS(a2),:);
                LeftTSS.Window(a,1:2) = x.Window(leftTSS(a2),1:2);
                LeftTSS.Enrichment(a,1) = x.Enrichment(leftTSS(a2),1);
                a = a + 1;
            end
            if isempty(a1) == 0 && isempty(a2) == 0
                LeftTSS.Window(a,1:2) = [x.Window(leftTSS(a1),1) x.Window(leftTSS(a2),2)];
                LeftTSS.Distance(a,:) = [mean(LeftTSS.Window(a,1:2)) i (mean(LeftTSS.Window(a,1:2)) - knownGene{i,start})];
                if j == 1 || a == 1
                    LeftTSS.Reads(a,1) = sum(x.Reads(leftTSS(a1:a2),1));
                    LeftTSS.Enrichment(a,1) = LeftTSS.Reads(a,1)/sum(PI(leftTSS(a1:a2),1));
                    a = a + 1;
                elseif j > 1 && (LeftTSS.Window(a,1) >= LeftTSS.Window(a-1,2)) 
                    LeftTSS.Reads(a,1) = sum(x.Reads(leftTSS(a1:a2),1));
                    LeftTSS.Enrichment(a,1) = LeftTSS.Reads(a,1)/sum(PI(leftTSS(a1:a2),1));
                    a = a + 1;
                else
                    LeftTSS.Window(a,:) = NaN;
                    LeftTSS.Distance(a,:) = NaN;
                end
            end
        end
        if a2 == length(distleftTSS)
            break
        end
        r = a2 + 1;
    end
    
    r = 1;
    for j = 1:length(rightTSS)
        distrightTSS = x.Distance(rightTSS,1) - x.Distance(rightTSS(r),1);
        if length(distrightTSS) == 1
            RightTSS.Distance(b,:) = x.Distance(rightTSS(r),[1 4 5]);
            RightTSS.Reads(b,1) = x.Reads(rightTSS(r),:);
            RightTSS.Window(b,1:2) = x.Window(rightTSS(r),1:2);
            RightTSS.Enrichment(b,1) = x.Enrichment(rightTSS(r),1);
            b = b + 1;
            break
        else
            a1 = find((abs(distrightTSS) <= window_thresh),1,'first');
            a2 = find((abs(distrightTSS) <= window_thresh),1,'last');
            if a1 == a2
                a1 = [];
                RightTSS.Distance(b,:) = x.Distance(rightTSS(a2),[1 4 5]);
                RightTSS.Reads(b,1) = x.Reads(rightTSS(a2),:);
                RightTSS.Window(b,1:2) = x.Window(rightTSS(a2),1:2);
                RightTSS.Enrichment(b,1) = x.Enrichment(rightTSS(a2),1);
                b = b + 1;
            end
            if isempty(a1) == 0 && isempty(a2) == 0
                RightTSS.Window(b,1:2) = [x.Window(rightTSS(a1),1) x.Window(rightTSS(a2),2)];
                RightTSS.Distance(b,:) = [mean(RightTSS.Window(b,1:2)) i (mean(RightTSS.Window(b,1:2)) - knownGene{i,start})];
                if j == 1 || b == 1
                    RightTSS.Reads(b,1) = sum(x.Reads(rightTSS(a1:a2),1));
                    RightTSS.Enrichment(b,1) = RightTSS.Reads(b,1)/sum(PI(rightTSS(a1:a2),1));
                    b = b + 1;
                elseif j > 1 && (RightTSS.Window(b,1) >= RightTSS.Window(b-1,2))
                    RightTSS.Reads(b,1) = sum(x.Reads(rightTSS(a1:a2),1));
                    RightTSS.Enrichment(b,1) = RightTSS.Reads(b,1)/sum(PI(rightTSS(a1:a2),1));
                    b = b + 1;
                else
                    RightTSS.Window(b,:) = NaN;
                    RightTSS.Distance(b,:) = NaN;
                end
            end
        end
        if a2 == length(distrightTSS)
            break
        end
        r = a2 + 1;
    end
    
    r = 1;
    for j = 1:length(leftTES)
        distleftTES = x.Distance(leftTES,6) - x.Distance(leftTES(r),6);
        if length(distleftTES) == 1
            LeftTES.Distance(c,:) = x.Distance(leftTES(r),6:8);
            LeftTES.Reads(c,1) = x.Reads(leftTES(r),:);
            LeftTES.Window(c,1:2) = x.Window(leftTES(r),1:2);
            LeftTES.Enrichment(c,1) = x.Enrichment(leftTES(r),1);
            c = c + 1;
            break
        else
            a1 = find((abs(distleftTES) <= window_thresh),1,'first');
            a2 = find((abs(distleftTES) <= window_thresh),1,'last');
            if a1 == a2
                a1 = [];
                LeftTES.Distance(c,:) = x.Distance(leftTES(a2),6:8);
                LeftTES.Reads(c,1) = x.Reads(leftTES(a2),:);
                LeftTES.Window(c,1:2) = x.Window(leftTES(a2),1:2);
                LeftTES.Enrichment(c,1) = x.Enrichment(leftTES(a2),1);
                c = c + 1;
            end
            if isempty(a1) == 0 && isempty(a2) == 0
                LeftTES.Window(c,1:2) = [x.Window(leftTES(a1),1) x.Window(leftTES(a2),2)];
                LeftTES.Distance(c,:) = [mean(LeftTES.Window(c,1:2)) i (mean(LeftTES.Window(c,1:2)) - knownGene{i,endit})];
                if j == 1 || c == 1
                    LeftTES.Reads(c,1) = sum(x.Reads(leftTES(a1:a2),1));
                    LeftTES.Enrichment(c,1) = LeftTES.Reads(c,1)/sum(PI(leftTES(a1:a2),1));
                    c = c + 1;
                elseif j > 1 && (LeftTES.Window(c,1) >= LeftTES.Window(c-1,2))
                    LeftTES.Reads(c,1) = sum(x.Reads(leftTES(a1:a2),1));
                    LeftTES.Enrichment(c,1) = LeftTES.Reads(c,1)/sum(PI(leftTES(a1:a2),1));
                    c = c + 1;
                else
                    LeftTES.Window(c,:) = NaN;
                    LeftTES.Distance(c,:) = NaN;
                end
            end
        end
        if a2 == length(distleftTES)
            break
        end
        r = a2 + 1;
    end
    
    r = 1;
    for j = 1:length(rightTES)
        distrightTES = x.Distance(rightTES,6) - x.Distance(rightTES(r),6);
        if length(distrightTES) == 1
            RightTES.Distance(d,:) = x.Distance(rightTES(r),[6 9 10]);
            RightTES.Reads(d,1) = x.Reads(rightTES(r),:);
            RightTES.Window(d,1:2) = x.Window(rightTES(r),1:2);
            RightTES.Enrichment(d,1) = x.Enrichment(rightTES(r),1);
            d = d + 1;
            break
        else
            a1 = find((abs(distrightTES) <= window_thresh),1,'first');
            a2 = find((abs(distrightTES) <= window_thresh),1,'last');
            if a1 == a2
                a1 = [];
                RightTES.Distance(d,:) = x.Distance(rightTES(a2),[6 9 10]);
                RightTES.Reads(d,1) = x.Reads(rightTES(a2),:);
                RightTES.Window(d,1:2) = x.Window(rightTES(a2),1:2);
                RightTES.Enrichment(d,1) = x.Enrichment(rightTES(a2),1);
                d = d + 1;
            end
            if isempty(a1) == 0 && isempty(a2) == 0
                RightTES.Window(d,1:2) = [x.Window(rightTES(a1),1) x.Window(rightTES(a2),2)];
                RightTES.Distance(d,:) = [mean(RightTES.Window(d,1:2)) i (mean(RightTES.Window(d,1:2)) - knownGene{i,endit})];
                if j == 1 || d == 1
                    RightTES.Reads(d,1) = sum(x.Reads(rightTES(a1:a2),1));
                    RightTES.Enrichment(d,1) = RightTES.Reads(d,1)/sum(PI(rightTES(a1:a2),1));
                    d = d + 1;
                elseif j > 1 && (RightTES.Window(d,1) >= RightTES.Window(d-1,2)) 
                    RightTES.Reads(d,1) = sum(x.Reads(rightTES(a1:a2),1));
                    RightTES.Enrichment(d,1) = RightTES.Reads(d,1)/sum(PI(rightTES(a1:a2),1));
                    d = d + 1;
                else
                    RightTES.Window(d,:) = NaN;
                    RightTES.Distance(d,:) = NaN;
                end
            end
        end
        if a2 == length(distrightTES)
            break
        end
        r = a2 + 1;
    end
end

A = (isnan(LeftTSS.Distance(:,1)) == 0);
LeftTSS.Distance = LeftTSS.Distance(A,:);
LeftTSS.Reads = LeftTSS.Reads(A);
LeftTSS.Window = LeftTSS.Window(A,:);
LeftTSS.Enrichment = LeftTSS.Enrichment(A);

A = (isnan(RightTSS.Distance(:,1)) == 0);
RightTSS.Distance = RightTSS.Distance(A,:);
RightTSS.Reads = RightTSS.Reads(A);
RightTSS.Window = RightTSS.Window(A,:);
RightTSS.Enrichment = RightTSS.Enrichment(A);

A = (isnan(LeftTES.Distance(:,1)) == 0);
LeftTES.Distance = LeftTES.Distance(A,:);
LeftTES.Reads = LeftTES.Reads(A);
LeftTES.Window = LeftTES.Window(A,:);
LeftTES.Enrichment = LeftTES.Enrichment(A);

A = (isnan(RightTES.Distance(:,1)) == 0);
RightTES.Distance = RightTES.Distance(A,:);
RightTES.Reads = RightTES.Reads(A);
RightTES.Window = RightTES.Window(A,:);
RightTES.Enrichment = RightTES.Enrichment(A);


New = struct('LeftTSS',LeftTSS,'RightTSS',RightTSS,'LeftTES',LeftTES,'RightTES',RightTES);
end
