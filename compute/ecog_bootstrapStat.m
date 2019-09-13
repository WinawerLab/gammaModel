function [zstat,upci_base] = ecog_bootstrapStat(data_bs,base_bs,option_test)

zstat = zeros(size(data_bs,1),1);

if isequal(option_test,'zstat')
    
    for kk = 1:size(data_bs,1) % number of stimuli
        % these are means of each bootstrap: 
        xmn = data_bs(kk,:);
        ymn = base_bs;

        % now we want to compute statistics
        mn1 = mean(xmn);
        mn2 = mean(ymn);
        se1 = std(xmn);  % this is already like a "standard error"
        se2 = std(ymn);
        zstat(kk) = (mn1-mn2)./sqrt(se1.^2+se2.^2); % take the standard errors and combine them
    end
    
elseif isequal(option_test,'bootstat')
    
	base_ci = prctile(base_bs,[16 84]);
    for kk = 1:size(data_bs,1) % number of stimuli
        data_ci = prctile(data_bs(kk,:),[16 84]);
        if data_ci(1)>base_ci(2)
            zstat(kk) = 1;
        end
    end
    upci_base = base_ci(2);
end
