n = 1000;
cod = NaN(1,n);
codranked = cod;

figure('Position',[0 0 600 300])
subplot(1,2,1),hold on
for ii = 1:1000
    x = [zeros(1,45) 1:5];
    data = x+randn(1,50)*.2;
    pred = x+randn(1,50)*.2;
    cod(ii) = 1-sum((data-pred).^2)/sum(data-mean(data).^2);
    
    [~, dataranks] = sort(data);
    [~, predranks] = sort(pred);

    codranked(ii) = 1 - sum((dataranks-predranks).^2)/sum((dataranks-mean(predranks)).^2);

    
    plot(pred,data,'k.')
end
xlabel('prediction'),ylabel('data')
title('simulate prediction and data')

subplot(1,2,2)
histogram(cod, -1:.01:1);
hold on
histogram(codranked, -1:.01:1);
legend('cod data', 'cod ranks','Location','northwest')
xlabel('cod (R^2)')
ylabel('number of simulations')

% Save  Figure
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['simulate_rankorder']))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
        ['simulate_rankorder']))