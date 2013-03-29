%Define the coverage of a transcription start site for a 1kbp window on a
%single chromosome

function x = coverage(chr,tss,data)
delta = 20;
TSS = floor(tss/25);

%Create a dynamic variable that changes with respect to the initial
%functions inputs

coordinates = sprintf('data.chip.chr%s',chr);
coordinatesF = sprintf('data.chipF.chr%s',chr);
coordinatesR = sprintf('data.chipR.chr%s',chr);
values = eval(coordinates);
valuesF = eval(coordinatesF);
valuesR = eval(coordinatesR);
if TSS <= delta
    x(1,1) = sum(values(1:TSS+delta));
    x(1,2) = sum(valuesF(1:TSS+delta)) + sum(valuesR(1:TSS+delta));
elseif TSS >= (length(values)-delta)
    x(1,1) = sum(values(TSS-delta:end));
    x(1,2) = sum(valuesF(TSS-delta:end)) + sum(valuesR(TSS+delta:end));
else
    x(1,1) = sum(values(TSS-delta:TSS+delta));
    x(1,2) = sum(valuesF(TSS-delta:TSS+delta)) + sum(valuesR(TSS-delta:TSS+delta));
end
end