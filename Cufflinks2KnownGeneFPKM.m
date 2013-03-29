clc

x = importdata('../../RNA-Seq/ESC_Wildtype/genes.fpkm_tracking');
load ../../RNA-Seq/known2Ref.mm9.mat
load knownGene.mm9.mat
ESC_RNA = x.data;
ESC_RNA_text = x.textdata(2:end,:);
clear x

greater = find(ESC_RNA(:,1) > 1);
ESC_RNA(isnan(ESC_RNA(:,1))|isempty(ESC_RNA(:,1)),1) = 0;
lam = poissfit(ESC_RNA(greater,1));
FPKM_Percentile = zeros(1,length(ESC_RNA));
FPKM_Percentile(1,greater) = 1 - poisscdf(ESC_RNA(greater,1),lam);

FPKM_ESC = cell(length(knownGene),5);
for i = 1:length(knownGene)
   x = find(strcmp(knownGene(i,1),known2Ref(:,1)));
   if ~isempty(x)
       Ref = known2Ref(x,2);
       y = find(strcmp(ESC_RNA_text(greater,1),Ref));
       FPKM = ESC_RNA(greater(y),1);
       FPKM_Per = FPKM_Percentile(1,greater(y));
       if ~isempty(FPKM)
           FPKM_ESC(i,:) = {knownGene{i,1} Ref knownGene{i,3} FPKM FPKM_Per};
       else
           FPKM_ESC(i,:) = {knownGene{i,1} Ref knownGene{i,3} NaN NaN};
       end
   else
       FPKM_ESC(i,:) = {knownGene{i,1} NaN knownGene{i,3} NaN NaN};
   end
   
end