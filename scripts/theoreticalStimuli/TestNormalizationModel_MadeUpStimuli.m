
% Pretend contrast energy for the following hypothetical stimuli:
% rows: there are 9 locations (3x3)
% columns: there are 4 orientations

stim=zeros(8,9,4);

% Number of orientations:
% grating
% plaid
% squiggles

% grating:
stim(1,:,:) = repmat([1 0 0 0],9,1);
% plaid:
stim(2,:,:) = repmat([.5 0 .5 0],9,1);
% random:
temp_s = repmat([.5 0 .5 0],9,1);
stim(3,:,:) = reshape(temp_s(randperm(length(temp_s(:)))),9,4);
clear temp_s

% Variations across space
% size: small grating
stim(4,:,:) = stim(1,:,:);
stim(4,[1:4 6:9],:) = 0;

% pattern
stim(5,:,:) = stim(2,:,:);
stim(5,[1:4 6:9],:) = 0;

% squiggles
stim(6,:,:) = stim(3,:,:);
stim(6,[1:4 6:9],:) = 0;

% grating surrounded by other orientation
stim(7,:,:) = repmat([0 1 0 0],9,1);
stim(7,5,:) = [1 0 0 0];

% grating to noise 
stim(8,:,:) = zeros(size(stim(1,:,:)));
for kk = 1:size(stim,2)
    stim(8,kk,randperm(4,1)) = 1;
end

% % plot all contrast energies for 4 different orientations:
% figure('Position',[0 0 300 600])
% for kk = 1:size(stim,1)
%     thisStim = reshape(squeeze(stim(kk,:,:)),3,3,4);
%     for rr = 1:4
%         subplot(size(stim,1),4,4*kk - 4 + rr)
%         imagesc(thisStim(:,:,rr),[0 1])
%         axis square
%         axis off
%     end
% end
%%
cm = copper(150);

figure('Position',[0 0 300 600])
for kk = 1:size(stim,1)
    thisStim = reshape(squeeze(stim(kk,:,:)),3,3,4);
    for rr = 1:4
        subplot(size(stim,1),4,4*kk - 4 + rr),hold on
        for ii = 1:size(thisStim(:,:,rr),1)
            for jj = 1:size(thisStim(:,:,rr),2)
                plot(ii,jj,'.','MarkerSize',50,...
                    'Color',cm(thisStim(ii,jj,rr)*100+50,:))
            end
        end
        xlim([.5 3.5]),ylim([0.5 3.5])
        axis square, box on
        set(gca,'XTick',[],'YTick',[])
        set(gca,'Color',[0 0 0])
    end
end

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/normalization_testStim/testStim_contrastE'])
% print('-depsc','-r300',['./figures/normalization_testStim/testStim_contrastE'])


%% First normalization model

n_use = [.5 1 2 6];
R_out = zeros(size(stim,1),length(n_use));
r_max = 1;
c50 = .2;

for kk = 1:size(stim,1) % images
    thisImage = squeeze(stim(kk,:,:));
    imFilt1Sum = sum(thisImage,1); % Sum across space 
    % imFilt1Sum = the dot product between simple cell contrast energy and weights

    for ii = 1:length(n_use) % different n's to try
        n = n_use(ii);
        
        % Normalization
        c_rms = sqrt(sum(imFilt1Sum.^2));
        R_out(kk,ii) = r_max * (sum(imFilt1Sum.^n))./...
            (c50.^n + c_rms.^n);
    end
end

figure,hold on
color_plot = {'r','k','b','g','y'};
for ii = 1:length(n_use) % different n's to try
    plot(1:3,R_out(1:3,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(4:6,R_out(4:6,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(7,R_out(7,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(R_out(:,ii),'.','Color',color_plot{ii},'MarkerSize',10)
end
xlim([0 size(stim,1)+1])

%% Why does this work?

n_use = [.5 1 2 4];
R_out = zeros(size(stim,1),length(n_use));
D_out = zeros(size(stim,1),length(n_use));
N_out = zeros(size(stim,1),length(n_use));
r_max = 1;
c50 = .2;

for kk = 1:size(stim,1) % images
    thisImage = squeeze(stim(kk,:,:));
%     imFilt1Sum = sum(thisImage,1); % Sum across space 
    imFilt1Sum = mean(thisImage,1); % Mean across space: do not take size in drive...
    % imFilt1Sum = the dot product between simple cell contrast energy and weights

    for ii = 1:length(n_use) % different n's to try
        n = n_use(ii);
        
        % Normalization
        c_rms = sqrt(sum(imFilt1Sum.^2));
        
        R_out(kk,ii) = r_max * (sum(imFilt1Sum.^n))./...
            (c50.^n + c_rms.^n);
        
        D_out(kk,ii) = (sum(imFilt1Sum.^n));
        N_out(kk,ii) = c_rms.^n;
    end
end

figure
subplot(3,1,1),hold on
color_plot = {'r','k','b','g','y'};
for ii = 1:length(n_use) % different n's to try
    plot(1:3,R_out(1:3,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(4:6,R_out(4:6,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(7,R_out(7,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(R_out(:,ii),'.','Color',color_plot{ii},'MarkerSize',10)
end
xlim([0 8])
title('Response')

subplot(3,1,2),hold on
color_plot = {'r','k','b','g','y'};
for ii = 1:length(n_use) % different n's to try
    plot(1:3,D_out(1:3,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(4:6,D_out(4:6,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(7,D_out(7,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(D_out(:,ii),'.','Color',color_plot{ii},'MarkerSize',10)
end
xlim([0 8])
title('Drive')

subplot(3,1,3),hold on
color_plot = {'r','k','b','g','y'};
for ii = 1:length(n_use) % different n's to try
    plot(1:3,N_out(1:3,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(4:6,N_out(4:6,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(7,N_out(7,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(N_out(:,ii),'.','Color',color_plot{ii},'MarkerSize',10)
end
xlim([0 8])
title('Normalizing term')



%% Center drive, all normalization pool
n_use = [.1 .5 1 2 3];
R_out = zeros(size(stim,1),length(n_use));
r_max = 1;
c50 = .2;

for kk = 1:size(stim,1) % images
    thisImage = squeeze(stim(kk,:,:));
    
    % Driving pool
    imFilt1Sum = sum(thisImage(5,:),1); % Sum across space
    % Normalization pool
    imFilt2Sum = sum(thisImage,1); % Sum across space
    
    for ii = 1:length(n_use) % different n's to try
        n = n_use(ii);
        
        % Normalization
        c_rms = sqrt(sum(imFilt2Sum.^2));
%         R_out(kk,ii) = r_max * (sum(imFilt1Sum.^n)) .* ...
%             ((c50.^n + c_rms.^n).^-1);
        R_out(kk,ii) = r_max * (sum(imFilt1Sum.^n));
    end
end

figure,hold on
color_plot = {'r','k','b','g','y'};
for ii = 1:length(n_use) % different n's to try
    plot(1:3,R_out(1:3,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(4:6,R_out(4:6,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(7,R_out(7,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(R_out(:,ii),'.','Color',color_plot{ii},'MarkerSize',10)
end

xlim([0 8])

%%
%% Center drive, normalization on/off, dependent on image homogeneity
n_use = [.5 1 2 4];
R_out = zeros(size(stim,1),length(n_use));
r_max = 1;
c50 = .2;

gatingTerm = zeros(size(stim,1),1);
beta_cntr = 1;

for kk = 1:size(stim,1) % images
    thisImage = squeeze(stim(kk,:,:));
    
    % Driving pool
    imFilt1Sum = sum(thisImage,1); % Sum across space 
    
    gatingTerm(kk) = sum(var(thisImage,[],1)); 

    for ii = 1:length(n_use) % different n's to try
        n = n_use(ii);
        
        % Normalization
        c_rms = sqrt(sum(imFilt1Sum.^2));
%         R_out(kk,ii) = r_max * (sum(imFilt1Sum.^n));
        
        if gatingTerm(kk)<.1 % yes surround suppression
            R_out(kk,ii) = r_max * (sum(imFilt1Sum.^n)) .* ...
                ((c50.^n + c_rms.^n).^-1);
        else % no surround suppression
            R_out(kk,ii) = r_max * (sum(imFilt1Sum.^n)) .* ...
                ((c50.^n + c_rms.^n).^-1);
        end
    end
end

%%
figure,hold on
color_plot = {'r','k','b','g','y'};
for ii = 1:length(n_use) % different n's to try
    plot(1:3,R_out(1:3,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(4:6,R_out(4:6,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(7,R_out(7,ii),'Color',color_plot{ii},'LineWidth',1)
    plot(R_out(:,ii),'.','Color',color_plot{ii},'MarkerSize',10)
end
xlim([0 8])
