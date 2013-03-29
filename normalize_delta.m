function normed = normalize_delta(yy)

%Normalize with respect to the background. In the even that the background
%is 0, then just save the enrichment value from the antibody

[m,n] = size(yy)
for i = 1:n-1
    normed(:,i) = yy(:,i) - yy(:,4);
end
