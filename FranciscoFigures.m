WithinWindow_FH_PI = [sum(abs(PI_significant.Distance) <= 1000) sum(abs(PI_significant.Distance) > 1000 & abs(PI_significant.Distance) <= 10000)];
WithinWindow_FH = [sum(abs(TBP_significant.Distance) <= 1000) sum(abs(TBP_significant.Distance) > 1000 & abs(TBP_significant.Distance) <= 10000)];
WithinWindow_HY = [sum(abs(TAF7L_significant.Distance) <= 1000) sum(abs(TAF7L_significant.Distance) > 1000 & abs(TAF7L_significant.Distance) <= 10000)];

figure(2)
subplot(4,3,1)
pie(WithinWindow_FH)
title('TBP')
subplot(4,3,2)
pie(WithinWindow_FH_PI)
title('PI')
subplot(4,3,3)
pie(WithinWindow_HY)
title('TAF7L')
legend('<= 1kb','1-10kb')
subplot(4,3,4)
boxplot(TBP1)
ylabel('1kb')
subplot(4,3,7)
boxplot(TBP10)
ylabel('1-10kb')
subplot(4,3,10)
bar(TBP_mean), hold on
errorbar(TBP_mean,TBP_std,'kx'), hold off
subplot(4,3,5)
boxplot(PI1)
subplot(4,3,8)
boxplot(PI10)
subplot(4,3,11)
bar(PI_mean), hold on
errorbar(PI_mean,PI_std,'kx'), hold off
subplot(4,3,6)
boxplot(TAF7L1)
subplot(4,3,9)
boxplot(TAF7L1)
subplot(4,3,12)
bar(TAF7L_mean), hold on
errorbar(TAF7L_mean,TAF7L_std,'kx'), hold off
