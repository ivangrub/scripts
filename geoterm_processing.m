es1 = (ESCZscore3(:,1) >= 3 & ESCZscore3(:,2) <= 1.5);
es2 = (ESCZscore3(:,1) <= 1.5 & ESCZscore3(:,2) >= 3);
es3 = (ESCZscore3(:,1) >= 3 & ESCZscore3(:,3) <= 1.5);
es4 = (ESCZscore3(:,1) <= 1.5 & ESCZscore3(:,3) >= 3);
es5 = (ESCZscore3(:,3) >= 3 & ESCZscore3(:,2) <= 1.5);
es6 = (ESCZscore3(:,3) <= 1.5 & ESCZscore3(:,2) >= 3);

tbptaf = unique(restgene(es1,1));
taftbp = unique(restgene(es2,1));
tbppol = unique(restgene(es3,1));
poltbp = unique(restgene(es4,1));
poltaf = unique(restgene(es5,1));
tafpol = unique(restgene(es6,1));