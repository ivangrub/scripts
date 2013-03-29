x = ESCZscore1(:,1) >= 3;
y = ESCZscore1(:,2) >= 3;
z = ESCZscore1(:,3) >= 3;

m = ESCZscore1(:,1) <= 0;
n = ESCZscore1(:,2) <= 0;
o = ESCZscore1(:,3) <= 0;

A = [x-m y-n z-o];

c1 = (A(:,1) == -1 & A(:,2) == -1 & A(:,3) == -1);
c2 = (A(:,1) == -1 & A(:,2) == -1 & A(:,3) == 1);
c3 = (A(:,1) == -1 & A(:,2) == 1 & A(:,3) == -1);
c4 = (A(:,1) == -1 & A(:,2) == 1 & A(:,3) == 1);
c5 = (A(:,1) == 1 & A(:,2) == -1 & A(:,3) == -1);
c6 = (A(:,1) == 1 & A(:,2) == -1 & A(:,3) == 1);
c7 = (A(:,1) == 1 & A(:,2) == 1 & A(:,3) == -1);
c8 = (A(:,1) == 1 & A(:,2) == 1 & A(:,3) == 1);
B= [sum(c1),sum(c2),sum(c3),sum(c4),sum(c5),sum(c6),sum(c7),sum(c8)];


expTBP3 = sum(x)/49842;
expTAF3 = sum(y)/49842;
expPol3 = sum(z)/49842;

expTBP0 = sum(m)/49842;
expTAF0 = sum(n)/49842;
expPol0 = sum(o)/49842;

expTBP = 1 - expTBP3-expTBP0;
expTAF = 1 - expTAF3-expTAF0;
expPol = 1 - expPol3-expPol0;

expB = [expTBP0*expTAF0*expPol0*49842 expTBP3*expTAF0*expPol0*49842 expTBP0*expTAF3*expPol0*49842 ...
    expTBP0*expTAF3*expPol3*49842 expTBP3*expTAF0*expPol0*49842 expTBP3*expTAF0*expPol3*49842 ...
    expTBP3*expTAF3*expPol0*49842 expTBP3*expTAF3*expPol3*49842];
    
C = [B' expB'];

