function normed = normalize_delta_PreImmune(yy,control)

%Normalize with respect to the background. In the even that the background
%is 0, then just save the enrichment value from the antibody

[~,n] = size(yy);
normed = zeros(control-1,n);

for i = 1:(control - 1)
   normed(i,:) = yy(i,:) - yy(control,:);
end
  