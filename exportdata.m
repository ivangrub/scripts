%Export the ChIP-Seq Data

function exportdata(knownGene,CpG,CellType,x)

str = date;
A = {knownGene{:,1};knownGene{:,2};knownGene{:,3};knownGene{:,4}}';
if strcmp(x,'ESC')
    headers = {'Accession Number' 'Gene' 'Description' 'Chromosome Number' ...
    'Forward/Reverse' 'TSS' 'TES' 'CpG Count' 'CpG Obs/Exp' 'TBP ChIP' ...
    'TBP ChIP F+R' 'TAF1 ChIP' 'TAF1 ChIP F+R' 'PolII_Ser5 ChIP' 'PolII_Ser5 ChIP F+R'...
    'PolII_N ChIP' 'PolII_N ChIP F+R' 'PI ChIP' 'PI ChIP F+R' 'MCPI ChIP' 'MCPI ChIP F+R'};
    xlswrite(sprintf('%s_CompiledData_%s',x,str),headers,'Compiled Data','A1');
    xlswrite(sprintf('%s_CompiledData_%s',x,str),A,'Compiled Data','A2');
    xlswrite(sprintf('%s_CompiledData_%s',x,str),knownGene(:,5:7),'Compiled Data','E2');
    xlswrite(sprintf('%s_CompiledData_%s',x,str),CpG,'Compiled Data','H2');
    xlswrite(sprintf('%s_CompiledData_%s',x,str),CellType,'Compiled Data','J2')
else
    headers = {'Accession Number' 'Gene' 'Description' 'Chromosome Number' ...
    'Forward/Reverse' 'TSS' 'TES' 'CpG Count' 'CpG Obs/Exp' 'TBP ChIP' ...
    'TBP ChIP F+R' 'PolII_Ser5 ChIP' 'PolII_Ser5 ChIP F+R' 'PI ChIP' 'PI ChIP F+R'...
    'Input ChIP' 'Input ChIP F+R'};
    xlswrite(sprintf('%s_CompiledData_%s',x,str),headers,'Compiled Data','A1');
    xlswrite(sprintf('%s_CompiledData_%s',x,str),A,'Compiled Data','A2');
    xlswrite(sprintf('%s_CompiledData_%s',x,str),knownGene(:,5:7),'Compiled Data','E2');
    xlswrite(sprintf('%s_CompiledData_%s',x,str),CpG,'Compiled Data','H2');
    xlswrite(sprintf('%s_CompiledData_%s',x,str),CellType,'Compiled Data','J2');
end

end
