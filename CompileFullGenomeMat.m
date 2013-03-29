%Compile a .mat file with all of the full genome variables

directory = {'Francisco'};
prefile = {'ESC'};
headers = {'TBP' 'PolII_Ser5' 'TAF1' 'PI'};
control = 4;

new_coverage_wholegenome;
chr_save;
Dist2TSS;
x = date;
save(sprintf('%s/%s_FullGenome_%s',directory,prefile,x))