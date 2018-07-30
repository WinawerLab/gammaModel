
% Pretend contrast energy for the following hypothetical stimuli:
% rows: there are 9 locations (3x3)
% columns: there are 4 orientations

stim=zeros(30,25,4);

% Contrast, Noise, Size
% As SOC:
% grating
% plaid
% squiggles

% grating, contrast:
stim(1,:,:) = repmat([1 0 0 0],25,1);
stim(2,:,:) = repmat([.8 0 0 0],25,1);
stim(3,:,:) = repmat([.6 0 0 0],25,1);
stim(4,:,:) = repmat([.4 0 0 0],25,1);
stim(5,:,:) = repmat([.2 0 0 0],25,1);

% plaid, contrast:
stim(6,:,:) = repmat([.5 0 .5 0],25,1);
stim(7,:,:) = repmat([.4 0 .4 0],25,1);
stim(8,:,:) = repmat([.3 0 .3 0],25,1);
stim(9,:,:) = repmat([.2 0 .2 0],25,1);
stim(10,:,:) = repmat([.1 0 .1 0],25,1);

% grating to noise, contrast
stim(11,:,:) = zeros(size(stim(1,:,:)));
for kk = 1:size(stim,2)
    stim(11,kk,randperm(4,1)) = 1;
end
stim(12,:,:) = .8*stim(11,:,:);
stim(13,:,:) = .6*stim(11,:,:);
stim(14,:,:) = .4*stim(11,:,:);
stim(15,:,:) = .2*stim(11,:,:);

% Size: smaller grating, plaid, noise
for kk = 1:15
    stim(15+kk,:,:) = stim(kk,:,:);
    stim(15+kk,[1:6 10 11 15 16 20 21:25],:) = 0;
end

%%
cm = copper(150);

figure('Position',[0 0 1200 200])
for kk = 1:size(stim,1)
    thisStim = reshape(squeeze(stim(kk,:,:)),5,5,4);
    for rr = 1:4
        subplot(4,size(stim,1),(rr-1)*30 + kk),hold on
        for ii = 1:size(thisStim(:,:,rr),1)
            for jj = 1:size(thisStim(:,:,rr),2)
                plot(ii,jj,'.','MarkerSize',15,...
                    'Color',cm(thisStim(ii,jj,rr)*100+50,:))
            end
        end
        xlim([.5 5.5]),ylim([0.5 5.5])
        axis square, box on
        set(gca,'XTick',[],'YTick',[])
        set(gca,'Color',[0 0 0])
    end
end

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/testStims_contrastE'])
% print('-depsc','-r300',['./figures/testStims_contrastE'])


%% Normalization model 1:

n_use = [.5 1 2];
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

figure('Position',[0 0 1250 400])
subplot(3,1,1),hold on
plot(R_out)
xlim([0 size(stim,1)+1])
title('Response')
set(gca,'XTick',[1:size(stim,1)]),grid on

subplot(3,1,2),hold on
plot(D_out)
xlim([0 size(stim,1)+1])
title('Drive')
set(gca,'XTick',[1:size(stim,1)]),grid on

subplot(3,1,3),hold on
plot(N_out)
xlim([0 size(stim,1)+1])
title('Normalizing term')
set(gca,'XTick',[1:size(stim,1)]),grid on

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/normalizationModel1'])
% print('-depsc','-r300',['./figures/normalizationModel1'])


%% Normalization model 2 (center/surround):

n_use = [.5 1 2];
R_out = zeros(size(stim,1),length(n_use));
D_out = zeros(size(stim,1),length(n_use));
N_out = zeros(size(stim,1),length(n_use));
r_max = 1;
c50 = .2;

for kk = 1:size(stim,1) % images
    thisImage = squeeze(stim(kk,[7:9 12:14 17:19],:));
%     imFilt1Sum = sum(thisImage,1); % Sum across space 
    imFilt1Sum = mean(thisImage,1); % Mean across space: the dot product between simple cell contrast energy and weights

    thisImage = squeeze(stim(kk,:,:));
    imFilt2Sum = mean(thisImage,1);

    for ii = 1:length(n_use) % different n's to try
        n = n_use(ii);
        
        % Normalization
        c_rms = sqrt(sum(imFilt2Sum.^2));
        
        R_out(kk,ii) = r_max * (sum(imFilt1Sum.^n))./...
            (c50.^n + c_rms.^n);
        
        D_out(kk,ii) = (sum(imFilt1Sum.^n));
        N_out(kk,ii) = c_rms.^n;
    end
end

figure('Position',[0 0 1250 400])
subplot(3,1,1),hold on
plot(R_out)
xlim([0 size(stim,1)+1])
title('Response')
set(gca,'XTick',[1:size(stim,1)]),grid on

subplot(3,1,2),hold on
plot(D_out)
xlim([0 size(stim,1)+1])
title('Drive')
set(gca,'XTick',[1:size(stim,1)]),grid on

subplot(3,1,3),hold on
plot(N_out)
xlim([0 size(stim,1)+1])
title('Normalizing term')
set(gca,'XTick',[1:size(stim,1)]),grid on

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/normalizationModel2'])
% print('-depsc','-r300',['./figures/normalizationModel2'])


