function remove_repeat_sequence(sequence,coordinates)

[m,~] = size(coordinates);

for i = 1:m
    sequence.Sequence(coordinates(i,1):coordinates(i,2)) = 'N';
end
sequence.Header = 'chr2L_Benjamin';

fastawrite('~/Desktop/chr2L_Benjamin.fa',sequence)