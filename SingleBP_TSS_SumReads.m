headers ={'TBP' 'TAF1'};
person = {'Lea' 'LW'};

load knownGene.mm9.mat
for i = 1:length(knownGene)
  if knownGene{i,5} == '+'
    index(:,i) = [cell2mat(knownGene(i,6))-400; ...
                  cell2mat(knownGene(i,6))+400]
  else
    index(:,i) = [cell2mat(knownGene(i,7))-400; ...
                  cell2mat(knownGene(i,7))+400]
  end
end

A = zeros(800,length(knownGene),length(headers));

for i = 1:length(headers)
  fid = fopen(sprintf('../ChIP-Seq Data/%s/%s_%s.mm9.bw',person{1},person{2}, ...
                    headers{i}))
  j = 0;
  while 1 
    j = j + 1;
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end 
    [~,strand,CHR,st,seq,~,~] = strread(tline,'%s%c%s%d%s%s%d', ...
                                        'delimiter','\t');
    st = st + 1;
    idx = st >= index(1,strcmp(knownGene(:,4),CHR{1})) & st < index(2,strcmp(knownGene(:,4),CHR{1}))
    if sum(idx) ~= 0
    	ind = find(idx);	
      	for k = 1:length(ind)
        	b = st - index(1,ind(k));
        	A(b:b+length(seq),ind(k),i) = A(b:b+length(seq)-1,ind(k),i) + 1;
      	end
    else
      	continue
    end
  end
  fclose(fid);clear fid
end

    
        