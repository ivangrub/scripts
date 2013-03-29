%Compare condensed windows from different analysis methods

Methods = {'Log' 'Poisson' 'Macs'};

for j = 1:length(Methods)
    for k = j+1:length(Methods)
        for i = 1:length(YY)
            left_left = YY(i,1) <= x(:,1);  
            left_left2 = YY(i,2) > x(:,1);
            right_left = YY(i,1) > x(:,1);
            right_left2 = YY(i,2) <= x(:,2);
            if left_left && left_left2
            elseif right_left && right_left2
            end
        end
    end
end