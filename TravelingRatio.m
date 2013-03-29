%Calculate the traveling ratio

fullcover = zeros(2,length(knownGene));

for k = 1:length(headers)
    x = eval(headers{k});
    chr = fieldnames(x.bp);
    for i = 1:find(strcmp(chr{end},knownGene(:,4)),1,'last')
        c = strcmp(knownGene(i,4),chr);
        if knownGene{i,5} == 1
            promoter = [cell2mat(knownGene(i,6))-30 cell2mat(knownGene(i,6))+300];
            rest = cell2mat(knownGene(i,7));
            
            p = x.bp.(chr{c})(1,:) >= promoter(1) & x.bp.(chr{c})(1,:)+read_length < promoter(2);
            pl = x.bp.(chr{c})(1,:) <= promoter(1) & x.bp.(chr{c})(1,:)+read_length < promoter(2);
            pl_ratio = (x.bp.(chr{c})(1,pl)+read_length - promoter(1))./read_length;
            pr = x.bp.(chr{c})(1,:) < promoter(2) & x.bp.(chr{c})(1,:)+read_length > promoter(2);
            pr_ratio = (promoter(2) - x.bp.(chr{c})(1,pr))./read_length;
            r = x.bp.(chr{c})(1,:) >= promoter(2) & x.bp.(chr{c})(1,:)+read_length <= rest;
            rl = x.bp.(chr{c})(1,:) <= promoter(2) & x.bp.(chr{c})(1,:)+read_length > promoter(2); 
            rl_ratio = (x.bp.(chr{c})(1,rl)+read_length - promoter(2))./read_length;
            rr = x.bp.(chr{c})(1,:) < rest & x.bp.(chr{c})(1,:)+read_length > rest;
            rr_ratio = (rest - x.bp.(chr{c})(1,rr))./read_length;
            
            
            fullcover(1:2,i) =[sum([p pl_ratio pr_ratio]);sum([r rr_ratio rl_ratio])];
        else
            promoter = [cell2mat(knownGene(i,7))+30 cell2mat(knownGene(i,7))-300];
            rest = cell2mat(knownGene(i,6));
            
            p = x.bp.(chr{c})(1,:)+read_length <= promoter(1) & x.bp.(chr{c})(1,:) > promoter(2);
            pr = x.bp.(chr{c})(1,:)+read_length >= promoter(1) & x.bp.(chr{c})(1,:) > promoter(2);
            pr_ratio = (promoter(1) - x.bp.(chr{c})(1,pr))./read_length;
            pl = x.bp.(chr{c})(1,:) < promoter(2) & x.bp.(chr{c})(1,:)+read_length > promoter(2);
            pl_ratio = (x.bp.(chr{c})(1,pl)+read_length - promoter(2))./read_length;
            r = x.bp.(chr{c})(1,:)+read_length <= promoter(2) & x.bp.(chr{c})(1,:) >= rest;
            rr = x.bp.(chr{c})(1,:) <= promoter(2) & x.bp.(chr{c})(1,:)+read_length > promoter(2);
            rr_ratio = (promoter(2) - x.bp.(chr{c})(1,rr))./read_length;
            rl = x.bp.(chr{c})(1,:) < rest & x.bp.(chr{c})(1,:)+read_length > rest;
            rl_ratio = (x.bp.(chr{c})(1,rl) + read_length - rest)./read_length;
            
            fullcover(1:2,i) =[sum([p pl_ratio pr_ratio]);sum([r rr_ratio rl_ratio])];
        end
    end
    ratio = fullcover(1,fullcover(1,:)>fullcover(2,:))./fullcover(2,fullcover(1,:)>fullcover(2,:));
    inv_ratio = -(fullcover(2,fullcover(2,:)>fullcover(1,:))./fullcover(1,fullcover(2,:)>fullcover(1,:)));
    ratio(ratio == Inf) = 0;
    ratio(isnan(ratio)) = fullcover(1,isnan(ratio));
    inv_ratio(inv_ratio == Inf) = 0;
    inv_ratio(isnan(inv_ratio)) = fullcover(2,isnan(inv_ratio));
    assignin('base',sprintf('%s_cover',headers{k}),fullcover);
    assignin('base',sprintf('%s_ratio',headers{k}),ratio);
    assignin('base',sprintf('%s_inv_ratio',headers{k}),inv_ratio);
end



