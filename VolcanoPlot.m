%Make a Volcano Plot of some data 

%X-axis is by definition the Log2(Fold Change)
%Y-axis is by definition -Log10(Pvalue)

control = find(strcmp('PI',headers));
reads_thresh = 0;
Pvalue = zeros(length(Reads),control-1);

A = normalize_fold(Reads,control,reads_thresh,'Lin');
Logged = log2(A);

for i = 1:control - 1
   [h,p] = ztest(Reads(:,[i control])); 
   Pvalue(:,i) = -log10(p);
end

plot(Logged,Pvalue)