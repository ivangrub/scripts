function normed = normalize_fold_PreImmune(yy,idx,control,thresh)

%Normalize with respect to the background. In the even that the background
%is 0, then just save the enrichment value from the antibody

[~,n] = size(yy);
normed = zeros(length(idx),n);

for i = 1:length(idx)
   normed(i,:) = yy(idx(i),:)./yy(control,:);
   for j = 1:n
       if (isnan(normed(i,j)) == 1) || (yy(idx(i),j) < thresh) || normed(i,j) == Inf
           normed(i,j) = 0;
       end
   end
end
  