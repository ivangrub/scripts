%Internal normalization

win = (FullGenome_Windows(1,2)-FullGenome_Windows(1,1))*25+24;
win_thresh = 500;
control = find(strcmp('PI',headers)); 

sig_coordinates = cell(1,3);

for k = 1:(control-1)
    r = (Reads(:,k) > 10);
    R = Reads(r,k);
    
    Win = [FullGenome_Windows(r,1)*25-24 FullGenome_Windows(r,2)*25];
    Distance = Dist(r,:);
    
    Prob = 1 - poisscdf(R,mean(R));
    P = Prob <= 1e-9;
    
    sig_coordinates(k,:) = {Win(P,1:2) Distance(P,:) R(P,1)};
    y.Window = sig_coordinates{k,1};
    y.Distance = sig_coordinates{k,2};
    y.Reads = sig_coordinates{k,3};
    
    B = length(y.Distance);
    Condensed = neighbors_fromwindows_internal(B,win_thresh,win,y,knownGene);
    
    assignin('base',sprintf('%s_Poiss_Sig',headers{k}),y)
    assignin('base',sprintf('%s_Poiss_Conden',headers{k}),Condensed)
end