PolIIPostEnrich = PolII_PostDiff_Coverage./PI_PostDiff_Coverage;
PolIIPostEnrich(isnan(PolIIPostEnrich)) = 0;
PolIIPostEnrich(PolIIPostEnrich==Inf) = PolII_PostDiff_Coverage(PolIIPostEnrich==Inf);

PolIIPreEnrich = PolII_PreDiff_Coverage./PI_PreDiff_Coverage;
PolIIPreEnrich(isnan(PolIIPreEnrich)) = 0;
PolIIPreEnrich(PolIIPreEnrich==Inf) = PolII_PreDiff_Coverage(PolIIPreEnrich==Inf);

MicroArrayEnrich = (1-cell2mat(CompiledMicroData(:,4)))-(1-cell2mat(CompiledMicroData(:,5)));
MicroArrayEnrich = MicroArrayEnrich';
% MicroArrayEnrich(isnan(MicroArrayEnrich)) = 0;
% MicroArrayEnrich(MicroArrayEnrich==Inf) = cell2mat(CompiledMicroData(MicroArrayEnrich==Inf,4));
idx = ~isnan(cell2mat(CompiledMicroData(:,5)));
PolIIDiff = PolIIPostEnrich - PolIIPreEnrich;

figure(1)
plot(PolIIPostEnrich(idx),1-cell2mat(CompiledMicroData(idx,4)),'r*',...
    PolIIPreEnrich(idx),1-cell2mat(CompiledMicroData(idx,5)),'b*')
xlabel('PolII Enrichment'),ylabel('Microarray Expression')

figure(2)
plot(PolIIDiff(idx),MicroArrayEnrich(idx),'r*')
xlabel('{\Delta}PolII Enrichment (Post vs Pre)')
ylabel('{\Delta}Microarray Expression (Post vs Pre)')
    