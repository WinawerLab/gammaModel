
% Correlate image contrast of patterns filtered with steerable pyramid
% filters with gamma/broadband power.
%
% Dora Hermes, 2017

addpath('~/Documents/git/ecogBasicCode/')
addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);
addpath(genpath('~/Documents/m-files/knkutils'));
addpath('/Users/dora/Documents/m-files/steer/')

%%
%% %%%%%% LOAD images %%%%%%%%
%%

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

% load images
load(fullfile(dataRootPath,'stimuli','task-soc_stimuli.mat'),'stimuli')

% save image energies
load(['./data/task-soc_stimuli_steerFilt01.mat'],'imEnergy')
load(['./data/task-soc_stimuli_steerFilt_high.mat'],'imEnergyHigh')
load(['./data/task-soc_stimuli_steerFilt_mid.mat'],'imEnergyMid')
load(['./data/task-soc_stimuli_steerFilt_low.mat'],'imEnergyLow')

% %%%%% Now to view an image, one can use:
% im_nr = 1;
% a = reshape(imEnergy(im_nr,:),dims(1),dims(2),8);
% figure,imagesc(a(:,:,4))

res = sqrt(size(imEnergy,2)/8);  % resolution of the pre-processed stimuli
imEnergyMean = zeros(size(imEnergy,1),res*res);
for kk = 1:size(imEnergy,1)
    thisImage = imEnergy(kk,:);
    thisImage = reshape(thisImage,res*res,8);
    imFilt1 = sum(thisImage,2); % Sum across orientations
    imEnergyMean(kk,:) = imFilt1;
end

% Switch order 2 and 3 (space and orientation) to match previous expected model input for NBF:
imEnergyResh = reshape(permute(reshape(imEnergy,size(imEnergy,1),[],8),[1 3 2]),size(imEnergy,1),res*res*8);
imEnergyReshHigh = reshape(permute(reshape(imEnergy,size(imEnergyHigh,1),[],8),[1 3 2]),size(imEnergy,1),res*res*8);
imEnergyReshMid = reshape(permute(reshape(imEnergy,size(imEnergyMid,1),[],8),[1 3 2]),size(imEnergy,1),res*res*8);
imEnergyReshLow = reshape(permute(reshape(imEnergy,size(imEnergyLow,1),[],8),[1 3 2]),size(imEnergy,1),res*res*8);

%% Load ECoG data

dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
%%%%% Pick a subject:
subjects = [19,23,24];
s = 3; subj = subjects(s);

%%%%% Pick an electrode:
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3
elec = 46;

% Choose an analysis type:
% analysisType = 'spectraRERP500';
% analysisType = 'spectra';
% analysisType = 'spectra500';
analysisType = 'spectra200';

% load the data:
dataFitName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
load(dataFitName,'resamp_parms')

% Broadband estimate, one value per image:
bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
ecog_bb = squeeze(median(resamp_parms(:,:,2),2))-bb_base;
ecog_bb_err=[squeeze(quantile(resamp_parms(:,:,2),.025,2)) ...
    squeeze(quantile(resamp_parms(:,:,2),.975,2))]'-bb_base;

ecog_g = squeeze(median(resamp_parms(:,:,3),2));
ecog_g_err=[squeeze(quantile(resamp_parms(:,:,3),.025,2)) ...
    squeeze(quantile(resamp_parms(:,:,3),.975,2))]';

ecog_a = squeeze(median(resamp_parms(:,:,6),2));
ecog_a_err=[squeeze(quantile(resamp_parms(:,:,6),.025,2)) ...
    squeeze(quantile(resamp_parms(:,:,6),.975,2))]';


%% Plot prf features

% Make a prf
res = sqrt(size(imEnergy,2)/8);  % resolution of the pre-processed stimuli
[~,xx,yy] = makegaussian2d(res,2,2,2,2);
pp = [ceil(res/2) ceil(res/2) ceil(res/10)];
% Define prf center Gaussian:
gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
    pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

im_nr_plot = [10];% 50 54 58];

for kk = im_nr_plot%:size(imsR,3) % plot features of this image
    figure('Position',[0 0 800 400])
    thisImage = imEnergy(kk,:);
    thisImage = reshape(thisImage,res*res,8);
    
    imFilt1 = thisImage .* repmat(gaufun1(pp),1,8);
%     imFilt2 = thisImage .* repmat(surroundfun(gaufun2(pp),gaufun1(pp)),1,8);
  
    ax1 = subplot(3,6,1);
    imagesc(stimuli(:,:,kk))  
    title(['Stimulus ' int2str(im_nr_plot)])
    axis off
    subplot(3,6,2); 
    imagesc(reshape(mean(thisImage,2),res,res)),hold on
    title('Contrast energy')
    subplot(3,6,3)
    imagesc(reshape(gaufun1(pp),res,res))
%     subplot(3,6,4)
%     imagesc(reshape(surroundfun(gaufun2(pp),gaufun1(pp)),res,res))

    for orientPlot = 1:8
        subplot(3,8,8+orientPlot)
        imagesc(reshape(imFilt1(:,orientPlot),res,res),[0 max(imFilt1(:))])
        axis square, axis off
    end
%     subplot(3,8,12)
%     title('Image within prf filtered by different orientations')
%     for orientPlot = 1:8
%         subplot(3,8,16+orientPlot)
%         imagesc(reshape(imFilt2(:,orientPlot),res,res),[0 max(imFilt2(:))])
%         axis square, axis off
%     end
%     subplot(3,8,20)
%     title('Image surrounding prf filtered by different orientations')
    
    colormap(ax1,'gray')
    
%     set(gcf,'PaperPositionMode','auto')
%     print('-dpng','-r300',['./figures/imageProperties/sub' subj '_el'  int2str(elec) '_image' int2str(kk)])
%     print('-depsc','-r300',['./figures/imageProperties/sub' subj '_el'  int2str(elec) '_image' int2str(kk)])

end

%%
%% Correlate contrast energy with ECoG for every pixel
%%

% space images to use for correlation:
img_use = [1:38];

corr_bb = zeros(1,size(imEnergyMean,2));
corr_g = zeros(1,size(imEnergyMean,2));
for kk = 1:size(imEnergyMean,2)
    if mod(kk,10000)==0
        disp(['pixel ' int2str(kk) ' of '  int2str(size(imEnergyMean,2))])
    end
    corr_bb(kk) = corr(imEnergyMean(img_use,kk),ecog_bb(img_use));
    corr_g(kk) = corr(imEnergyMean(img_use,kk),ecog_g(img_use));
end

figure('Position',[0 0 600 300])
subplot(1,2,1)
imagesc(reshape(fisherz(corr_bb),res,res),[-3 3]),hold on
axis square, axis off
colorbar
title('fisherz(r): CE and bb')

subplot(1,2,2)
imagesc(reshape(fisherz(corr_g),res,res),[-3 3]),hold on
axis square, axis off
colorbar
title('fisherz(r): CE and gamma')

colormap jet

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['./figures/corr_CE_ECoG/sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_corrEcogCe'])
print('-dpng','-r300',['./figures/corr_CE_ECoG/sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_corrEcogCe'])


%%
%% Normalization test 1

% S3 El46: 371.9099  375.1703    0.5229    0.2630    0.8879 

pp = [ceil(res/2) ceil(res/2) ceil(res/10)];
% pp = [100 100 .3];
pp = [371.9099/8  375.1703/8    0.5229];

[~,xx,yy] = makegaussian2d(100,2,2,2,2);
gaufun2 = @(pp) vflatten(makegaussian2d(100,pp(1),pp(2),pp(3),...
    pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

% im_nr_plot = [10 39 50 54 58];
im_nr_plot = [1:86];

n_use = [.5 1 2 4];
R_out = zeros(length(im_nr_plot),length(n_use));
D_out = zeros(length(im_nr_plot),length(n_use));
N_out = zeros(length(im_nr_plot),length(n_use));
r_max = 1;
c50 = .5;

for kk = 1:length(im_nr_plot) % images
%     thisImage = imEnergy(kk,:);
%     thisImage = reshape(thisImage,res*res,8);
    
    thisImageIn = imEnergyHigh(kk,:);
    thisImageIn = reshape(thisImageIn,res*res,8);
    thisImage = zeros(100*100,8);
    for ii = 1:8
        a = reshape(thisImageIn(:,ii),res,res);
        im = processmulti(@imresize,a,[100 100]);
        thisImage(:,ii) = im(:);
    end

    imFilt1 = thisImage .* repmat(gaufun2(pp),1,8);
    imFilt1Sum = sum(imFilt1,1); % Sum across space, 
    % imFilt1Sum = the dot product between simple cell contrast energy and weights

    for ii = 1:length(n_use) % different n's to try
        n = n_use(ii);
        
        % Normalization
        dd = (sum(imFilt1Sum.^n));
        nn = sqrt(sum(imFilt1Sum.^2)).^n;
        R_out(kk,ii) = r_max * dd .*...
            ((c50.^n + nn).^-1);
        D_out(kk,ii) = dd;
        N_out(kk,ii) = nn;
    end
end
clear dd nn

figure('Position',[0 0 1000 300])
for nn = 1:length(n_use)
    stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];

    subplot(length(n_use),3,3*nn-2),hold on
    plot(R_out(:,nn),'k','LineWidth',2)
    title(['n = ' num2str(n_use(nn),3) '; r_m_a_x = ' num2str(r_max,2) '; c_5_0 = ' num2str(c50)])
    for k=1:length(stim_change) % plot stimulus cutoffs
        plot([stim_change(k) stim_change(k)],[0 max(R_out(:,nn))],'Color',[.5 .5 .5],'LineWidth',1)
    end
    set(gca,'XTick',stim_change,'XTickLabel',{},'Color',[1 1 1],'Box','off')
    axis tight
    
    subplot(length(n_use),3,3*nn-1),hold on
    plot(D_out(:,nn),'k','LineWidth',2)
    for k=1:length(stim_change)% plot stimulus cutoffs
        plot([stim_change(k) stim_change(k)],[0 max(D_out(:,nn))],'Color',[.5 .5 .5],'LineWidth',1)
    end
    set(gca,'XTick',stim_change,'XTickLabel',{},'Color',[1 1 1],'Box','off')
    axis tight

    subplot(length(n_use),3,3*nn),hold on
    plot(N_out(:,nn),'k','LineWidth',2)
    % plot stimulus cutoffs
    for k=1:length(stim_change)
        plot([stim_change(k) stim_change(k)],[0 max(N_out(:,nn))],'Color',[.5 .5 .5],'LineWidth',1)
    end
    set(gca,'XTick',stim_change,'XTickLabel',{},'Color',[1 1 1],'Box','off')
    axis tight
end

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/normalization/model_patternCenter'])
% print('-depsc','-r300',['./figures/normalization/model_patternCenter'])


%%
%% Normalization test 2: pools small and large with surround on for homogeneous images

pp = [ceil(res/2) ceil(res/2) 80];

% im_nr_plot = [10 39 50 54 58];
im_nr_plot = [1:86];

n_use = [.5 1 2 4];
R_out = zeros(length(im_nr_plot),length(n_use));
D_out = zeros(length(im_nr_plot),length(n_use));
N_out = zeros(length(im_nr_plot),length(n_use));
r_max = 1;
c50 = .5;

pp1 = gaufun1(pp); %prf
pp2 = gaufun1(pp.* [1 1 2]); % prf 2x size
pp3 = gaufun1(pp.* [1 1 3]); % prf 3x size
pp4 = gaufun1(pp.* [1 1 4]); % prf 4x size

prf_center = repmat(pp1,1,8);

% Make a surround prf:
prf_surround = pp2 - .5*pp1; % (2*prf - 1*prf)
prf_surround(prf_surround<0)=0; % remove values < 0
prf_surround = prf_surround*(1./sum(prf_surround)); % set surface are to 1
prf_surround = repmat(prf_surround,1,8);

% Get a measure for periodicity in prf and surround
out_fourier = getFourierImageVals(double(stimuli),pp,res,xx,yy);
gatingTerm = out_fourier.rptPrf; % if larger than .5, lots of periodicity

for kk = 1:length(im_nr_plot) % images
    disp(['im nr' int2str(kk)])
    thisImage = imEnergy(kk,:);
    thisImage = reshape(thisImage,res*res,8);
    
    % Center pool
    imFilt1 = thisImage .* prf_center; % imFilt1Sum = the dot product between simple cell contrast energy and weights
    imFilt1Sum = sum(imFilt1,1); % Sum across space
    % Surround pool
    imFilt2 = thisImage .* prf_surround;
    imFilt2Sum = sum(imFilt2,1); % Sum across space
    
    for ii = 1:length(n_use) % different n's to try
        n = n_use(ii);
        
        if gatingTerm(kk) > 0.5 % image periodic
            nn = sqrt(sum(imFilt1Sum.^2 + imFilt2Sum.^2)).^n;
        else % image non-periodic
            nn = sqrt(sum(imFilt1Sum.^2)).^n;
        end
        
        % Normalization
        dd = (sum(imFilt1Sum.^n));
        
        R_out(kk,ii) = r_max * dd .*...
            ((c50.^n + nn).^-1);
        D_out(kk,ii) = dd;
        N_out(kk,ii) = nn;
    end
end
clear dd nn

figure('Position',[0 0 1000 300])
for nn = 1:length(n_use)
    stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];

    subplot(length(n_use),3,3*nn-2),hold on
    plot(R_out(:,nn),'k','LineWidth',2)
    title(['n = ' num2str(n_use(nn),3) '; r_m_a_x = ' num2str(r_max,2) '; c_5_0 = ' num2str(c50)])
    for k=1:length(stim_change) % plot stimulus cutoffs
        plot([stim_change(k) stim_change(k)],[0 max(R_out(:,nn))],'Color',[.5 .5 .5],'LineWidth',1)
    end
    set(gca,'XTick',stim_change,'XTickLabel',{},'Color',[1 1 1],'Box','off')
    axis tight
    
    subplot(length(n_use),3,3*nn-1),hold on
    plot(D_out(:,nn),'k','LineWidth',2)
    for k=1:length(stim_change)% plot stimulus cutoffs
        plot([stim_change(k) stim_change(k)],[0 max(D_out(:,nn))],'Color',[.5 .5 .5],'LineWidth',1)
    end
    set(gca,'XTick',stim_change,'XTickLabel',{},'Color',[1 1 1],'Box','off')
    axis tight

    subplot(length(n_use),3,3*nn),hold on
    plot(N_out(:,nn),'k','LineWidth',2)
    % plot stimulus cutoffs
    for k=1:length(stim_change)
        plot([stim_change(k) stim_change(k)],[0 max(N_out(:,nn))],'Color',[.5 .5 .5],'LineWidth',1)
    end
    set(gca,'XTick',stim_change,'XTickLabel',{},'Color',[1 1 1],'Box','off')
    axis tight
end

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/normalization/model_patternCenter'])
% print('-depsc','-r300',['./figures/normalization/model_patternCenter'])

%% Plot image features

pp = [ceil(res/2) ceil(res/2) ceil(res/10)];

% S3 El46: 371.9099  375.1703    0.5229    0.2630    0.8879 
pp = [371.9099  375.1703   5];

res = sqrt(size(imEnergy,2)/8);  % resolution of the pre-processed stimuli
[~,xx,yy] = makegaussian2d(res,2,2,2,2);

% Define prf center Gaussian:
gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
    pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

out_raw = getRawImageVals(double(stimuli),pp,res,xx,yy);
out_norm = getDivNormImageVals(imEnergy,pp,res,xx,yy);
out_ce = getCEImageVals(imEnergy,pp,res,xx,yy);
out_fourier = getFourierImageVals(double(stimuli),pp,res,xx,yy);

%% Figure image feature:

stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];

nr_plots = 5;

figure('Position',[0 0 800 700])
subplot(nr_plots,1,1),hold on
plot(out_ce.prf1CE,'r','LineWidth',2)
title('CE in prf')
ylim([0 1.1])

subplot(nr_plots,1,2),hold on
plot(out_ce.prfVarCE,'r','LineWidth',2)
title('Var in CE in prf')
% ylim([0 .5])

subplot(nr_plots,1,3),hold on
plot(out_norm.NBF,'r','LineWidth',2)
title('NBF')
ylim([0 .3])

subplot(nr_plots,1,4),hold on
plot(out_fourier.rptPrf,'r','LineWidth',2)
title('periodicity in Prf')
% ylim([0 1.5])

subplot(nr_plots,1,5),hold on
plot(out_norm.NormR,'LineWidth',2)
title('Divisive Normalization Models')
% ylim([0 2.5])

% subplot(nr_plots,1,6),hold on
% plot(out_raw.varIM,'r','LineWidth',2)
% title('contrast variance in image pixels')
% ylim([0 .5])
%%
for kk = 1:nr_plots
    subplot(nr_plots,1,kk)
    for k=1:length(stim_change)% plot stimulus cutoffs
        plot([stim_change(k) stim_change(k)],[0 1],'Color',[.5 .5 .5],'LineWidth',1)
    end
    set(gca,'XTick',stim_change,'XTickLabel',{},'Color',[1 1 1],'Box','off')
end

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/corr_CE_ECoG/model_patternCenter_severalValues'])
% print('-depsc','-r300',['./figures/corr_CE_ECoG/model_patternCenter_severalValues'])
