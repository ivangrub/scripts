left1 = abs(Pol2_Ser5_significant.Distance(:,1)) <= 1000;
right1 = abs(Pol2_Ser5_significant.Distance(:,2)) <= 1000;
left10 = (abs(Pol2_Ser5_significant.Distance(:,1)) > 1000 & abs(Pol2_Ser5_significant.Distance(:,1)) <= 10000);
right10 = (abs(Pol2_Ser5_significant.Distance(:,2)) > 1000 & abs(Pol2_Ser5_significant.Distance(:,2)) <= 10000);

left_1 = find(left1);
right_1 = find(right1);
left_10 = find(left10);
right_10 = find(right10);

Pol2_Ser51_left = 2.^(Pol2_Ser5_significant.Enrichment(left1));
Pol2_Ser51_right = 2.^(Pol2_Ser5_significant.Enrichment(right1));
Pol2_Ser510_left = 2.^(Pol2_Ser5_significant.Enrichment(left10));
Pol2_Ser510_right = 2.^(Pol2_Ser5_significant.Enrichment(right10));

Pol2_Ser5_mean(1:2,1) = [mean(Pol2_Ser51_left); mean(Pol2_Ser51_right)];
Pol2_Ser5_mean(1:2,2) = [mean(Pol2_Ser510_left); mean(Pol2_Ser510_right)];
Pol2_Ser5_std(1:2,1) = [std(Pol2_Ser51_left) std(Pol2_Ser51_right)];
Pol2_Ser5_std(1:2,2) = [std(Pol2_Ser510_left) std(Pol2_Ser510_right)];

Pol2_Ser5gene = {Pol2_Ser5_significant.Genes(left_1,1:4) Pol2_Ser5_significant.Genes(right_1,5:8) ...
    Pol2_Ser5_significant.Genes(left_10,1:4) Pol2_Ser5_significant.Genes(right_10,5:8)};

clear left1 right1 right10 left10