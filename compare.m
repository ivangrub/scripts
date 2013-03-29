%Compare ChIP-Seq data between 2 different datasets

categories = {'uniquegene' 'cpggene'  'restgene'};

for i = 1:3
    x = eval(sprintf('ESCZscore%d',i));
    y = eval(sprintf('MNZscore%d',i));
    
    %Find the TBP and PolII sites that are in ESC and not in the MN
    TBP1 = (x(:,1) >=3 & y(:,1) <= 3);  
    PolII1 = (x(:,3) >= 3 & y(:,2) <= 3);

    %Find the TBP and PolII sites that are in MN and not in the ESC
    TBP2 = (x(:,1) <=3 & y(:,1) >= 3); 
    PolII2 = (x(:,3) <= 3 & y(:,2) >= 3);
      
    j1 = find(TBP1);
    j2 = find(PolII1);
    j3 = find(TBP2);
    j4 = find(PolII2);
    
    %Get appropriate genes for each subcategory
    k = mod(i,3);
    if k == 0
        k = 3;
    end
    group = eval(categories{k});
    
    %Save into their specific matrices
    assignin('base',sprintf('InESC_NoMN_TBP_case%d',k),group(j1,1:2));
    assignin('base',sprintf('InESC_NoMN_PolII_case%d',k),group(j2,1:2));
    assignin('base',sprintf('InMN_NoESC_TBP_case%d',k),group(j3,1:2));   
    assignin('base',sprintf('InMN_NoESC_PolII_case%d',k),group(j4,1:2));
end