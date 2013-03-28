x = {'Strand_median' 'Strand_max' 'Strand_min' 'Median_dist' 'Max_dist' 'Min_dist'};
vartable = zeros([3,length(x)]);

for i = 1:length(x)
    [coeffs,score,latent] = pca(eval(sprintf('%s',x{i})));
    if strfind(x{i},'Strand') == 1
        subplot(2,3,i)
        vartable(1:3,i) = cumsum(latent)./sum(latent);
        biplot(coeffs(:,1:3),'score',score(:,1:3),'VarLabels',{'+ Expression','- Expression','Distance'})
        title(sprintf('%s',x{i}))
    else
        subplot(2,3,i)
        vartable(1:2,i) = cumsum(latent)./sum(latent);
        biplot(coeffs(:,1:2),'score',score(:,1:2),'VarLabels',{'Expression','Distance'})
        title(sprintf('%s',x{i}))
    end
end