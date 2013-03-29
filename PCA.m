function [PC chip_standard coeff latent cumvar] = PCA(chip_wreads)

%Runs the PCA for the Chip-Seq data

[m,n] = size(chip_wreads);

%Standardize the data using the zscore function keeping in mind that each
%column is independent of the others
for i = 1:n
    chip_standard(:,i) = zscore(chip_wreads(:,i));
end

%Calculate the variance and coefficients of the principle components (PC)
[coeff,PC,latent,tsquare] = princomp(chip_standard);

%Culmulative variance
cumvar = cumsum(latent)./sum(latent);

end