Case1(:,1) = ESCZscore1(:,1) >= 3;
Case1(:,2) = ESCZscore1(:,2) >= 3;
Case1(:,3) = ESCZscore1(:,3) >= 3;
Case2(:,1) = ESCZscore2(:,1) >= 3;
Case2(:,2) = ESCZscore2(:,2) >= 3;
Case2(:,3) = ESCZscore2(:,3) >= 3;
Case3(:,1) = ESCZscore3(:,1) >= 3;
Case3(:,2) = ESCZscore3(:,2) >= 3;
Case3(:,3) = ESCZscore3(:,3) >= 3;

Case1_TBPTAF1 = (Case1(:,1) == 1 & Case1(:,2) == 1);
Case1_TBPPol2 = (Case1(:,1) == 1 & Case1(:,3) == 1);
Case1_TAF1Pol2 = (Case1(:,2) == 1 & Case1(:,3) == 1);
Case1_TBPTAF1Pol2 = (Case1(:,1) == 1 & Case1(:,2) == 1 & Case1(:,3) == 1);
Case2_TBPTAF1 = (Case2(:,1) == 1 & Case2(:,2) == 1);
Case2_TBPPol2 = (Case2(:,1) == 1 & Case2(:,3) == 1);
Case2_TAF1Pol2 = (Case2(:,2) == 1 & Case2(:,3) == 1);
Case2_TBPTAF1Pol2 = (Case2(:,1) == 1 & Case2(:,2) == 1 & Case2(:,3) == 1);
Case3_TBPTAF1 = (Case3(:,1) == 1 & Case3(:,2) == 1);
Case3_TBPPol2 = (Case3(:,1) == 1 & Case3(:,3) == 1);
Case3_TAF1Pol2 = (Case3(:,2) == 1 & Case3(:,3) == 1);
Case3_TBPTAF1Pol2 = (Case3(:,1) == 1 & Case3(:,2) == 1 & Case3(:,3) == 1);


