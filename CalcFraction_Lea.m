headers = {'H4Ac' 'H3'};
cage_data = {'Cage_MinData_MacroSpecific'};

fid = fopen('KnownGene_wTSS.txt','r');

kg = textscan(fid,'%s%s%d%s%s%s%d%d','delimiter','\t');
fclose(fid);clear fid


gene = kg{4};
[unique_gene,index] = unique(gene);
chr = kg{5}(index);
tss = kg{7};
tes = kg{8};
strand = kg{6}(index);

frac_info = zeros(3*(length(headers)*length(cage_data)),length(unique_gene));

INDEX = 1;
for j = 1:length(cage_data)
    y = eval(cage_data{j});
    for i = 1:length(headers)
        x = eval(headers{i});
        for k = 1:length(unique_gene)
            ind = find(strcmp(unique_gene(k),gene));
            if ~isempty(strfind(chr{k},'random'))
                continue
            end
            cage = y.bp.(chr{k});
            if strcmp(strand{k},'+')
                tss1 = min(tss(ind))-250;
                tss2 = max(tss(ind))+250;
            else
                tss1 = min(tes(ind))-250;
                tss2 = max(tes(ind))+250;
            end
            vec = zeros(1,length(tss1:tss2));
            vec_cage = zeros(1,length(tss1:tss2));
            inside = find(x.bp.(chr{k}) >= tss1 & x.bp.(chr{k}) <= tss2);
            inside_cage = find(cage >= tss1 & cage < tss2);
            for n = 1:length(inside)
                coord = x.bp.(chr{k})(inside(n)) - tss1 + 1;
                vec(coord) = vec(coord) + 1;
            end
            for n = 1:length(inside_cage)
                coord = cage(inside_cage(n)) - tss1 + 1;
                vec_cage(coord) = vec_cage(coord) + 1;
            end
            total = sum(vec);
            total_cage = sum(vec_cage);
           % bins = round(length(vec)/10);
            ind = 1:25:length(vec);
            freq = zeros(1,length(ind));
            freq_cage = zeros(1,length(ind));
            
            for z = 1:length(ind)-1
                freq(z) = sum(vec(round(ind(z)):round(ind(z+1))));
            end
            
            for z = 1:length(ind)-1
                freq_cage(z) = sum(vec_cage(round(ind(z)):round(ind(z+1))));
            end
            frac = freq/total;
            frac_cage = freq_cage/total_cage;
            
            if sum(frac_cage >= 0.8) > 0
                type = 4;
            elseif sum(frac_cage >= 0.2 & frac_cage < 0.8) >= 2
                type = 3;
            elseif sum(frac_cage >= 0.2 & frac_cage < 0.8) == 1
                type = 2;
            elseif sum(frac_cage >= 0.2) == 0
                type = 1;
            end
            [yy,zz] = max(frac_cage);
            if zz <= 10
                frac_info(INDEX:INDEX+2,k) = [yy;max(frac(1:zz+10));type];
            elseif zz >= length(frac)-10
                frac_info(INDEX:INDEX+2,k) = [yy;max(frac(zz-10:end));type];
            else
                frac_info(INDEX:INDEX+2,k) = [yy;max(frac(zz-10:zz+10));type];
            end
        end
        INDEX = INDEX + 3;
    end
end
[m,~] =size(frac_info); 
k = 1;
colors = {'r.' 'g.' 'b.' 'k.'};
for i = 1:3:m
    figure(k)
    for j = 1:4
        t = frac_info(i+2,:) == j & ~isnan(frac_info(i,:)) & ~isnan(frac_info(i+1,:));
        plot(frac_info(i,t),frac_info(i+1,t),colors{j})
        hold on
    end
    hold off
    xlabel('Cage Fraction')
    ylabel(headers{k})
    k = k + 1;
end
