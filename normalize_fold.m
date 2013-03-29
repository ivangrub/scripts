function normed = normalize_fold(yy,control,thresh,scale)

%Normalize with respect to the background. In the even that the background
%is 0, then just save the enrichment value from the antibody

[m,n] = size(yy);
normed = zeros(m,control-1);
a = strcmp(scale,'Log');
if a == 1
    x = log2(thresh);
else
    x = thresh;
end

for i = 1:(control - 1)
   normed(:,i) = yy(:,i)./yy(:,control);
   for j = 1:m
       if (normed(j,i) == Inf) || (isnan(normed(j,i)) == 1) || (yy(j,control) <= x)
           normed(j,i) = 0;
       end
   end
end
  