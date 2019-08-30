function zstat = ecog_bootstrapStat(data_bs,base_bs)

zstat = zeros(size(data_bs,1),1);

for kk = 1:size(data_bs,1)
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