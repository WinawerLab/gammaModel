%
% Correlate image contrast of patterns filtered with gabor filters with
% gamma/broadband power.
%
% Dora Hermes, 2017

clear all
% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath('/Users/dora/Documents/git/ecogBasicCode/')
addpath(genpath('/Users/dora/Documents/m-files/knkutils'));

%% Load original images
origIm = load(fullfile(dataDir,'soc_bids','stimuli','task-soc_stimuli.mat'),'stimuli');

%% Load preprocessed images and divisive normalization:

load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

% compute the square root of the sum of the squares of the outputs of
% quadrature-phase filter pairs (this is the standard complex-cell energy model).
% after this step, stimulus is images x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2));

% % compute the population term in the divisive-normalization equation.
% % this term is simply the average across the complex-cell outputs
% % at each position (averaging across orientation).
% stimulusPOP = blob(stimulus,2,8)/8;
% 
% % repeat the population term for each of the orientations
% stimulusPOP = upsamplematrix(stimulusPOP,8,2,[],'nearest');
% 
% % apply divisive normalization to the complex-cell outputs.  there are two parameters
% % that influence this operation: an exponent term (r) and a semi-saturation term (s).
% % the parameter values specified here were determined through a separate fitting
% % procedure (see paper for details).  for the purposes of this script, we will
% % simply hard-code the parameter values here and not worry about attempting to fit
% % the parameters.
% r = 1;
% s = 0.5;
% stimulus = stimulus.^r ./ (s.^r + stimulusPOP.^r);
% clear stimulusPOP;

% sum across orientation.  after this step, stimulus is images x positions.
imEnergyMean = blob(stimulus,2,8);

%% Display NBF model results for one electrode

%%%%% Pick a subject:
subjects = [19,23,24];
s = 3; subj = subjects(s);

% %%%% Pick an electrode:
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3

analysisType = 'spectra200';
modelType = 'NBFsimple2';

elec = 45;
res = sqrt(size(stimulus,2)/8);

% Load NBF model fit:
load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveNBF',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'xys_pix','seed_params','prf_sizes',...
        'cross_NBFparams','cross_NBFestimate')

% Load ECoG data:
dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
    'ses-01','derivatives','ieeg',...
    ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
load(dataFitName)

% Calculate ECoG gamma power percent signal change:
bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
ecog_bb = mean(100*(10.^(resamp_parms(:,:,2)-bb_base)-1),2);
ecog_bb_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,2),.16,2)) ...
    squeeze(quantile(resamp_parms(:,:,2),.84,2))]'-bb_base)-1);
ecog_g = mean(100*(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1),2);
ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
    squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);
ylims = [min(ecog_g_err(:)) max([ecog_g_err(:)])];


res = sqrt(size(imEnergyMean,2));  % resolution of the pre-processed stimuli
[~,xx,yy] = makegaussian2d(res,2,2,2,2);

% Define a Gaussian centered at the prf:
gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
    pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

im_nrs = [10 74:78];
% im_nrs = [69:73];

% Plot filtered (contrast) images with prf

figure('Position',[0 0 800 600])

max_prf_images = 0.05;

for kk = 1:length(im_nrs)
    % images + prf
    subplot(3,length(im_nrs),kk)
    imagesc(double(origIm.stimuli(:,:,im_nrs(kk))),[0 max(origIm.stimuli(:))]),hold on
    % Move pRF in 135x135 to original size of 800x800 
    %   half of the screen is 60 in filtered and 400 in original
    orig_x = 400 + (seed_params(1)-res/2) * 400./60;
    orig_y = 400 + (seed_params(2)-res/2) * 400./60;
    orig_sigma = seed_params(3) * 400./60;
    axis image

    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*orig_sigma);
    plot(c.x + orig_y, c.y + orig_x, 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
    [c.x, c.y] = pol2cart(c.th, 2*ones(1,numPoints)*orig_sigma);
    plot(c.x + orig_y, c.y + orig_x, 'r:') % this is just reversed because plot and imagesc are opposite, checked this with contour
    

    xlim([orig_y-3*orig_sigma orig_y+3*orig_sigma])
    ylim([orig_x-3*orig_sigma orig_x+3*orig_sigma])
    title(['im ' int2str(im_nrs(kk))])
    
    
    % contrast images + prf
    subplot(3,length(im_nrs),length(im_nrs)+kk)
    imagesc(reshape(imEnergyMean(im_nrs(kk),:),res,res),[0 2]);
    axis image
    colormap gray
    hold on
    
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*seed_params(3));
    plot(c.x + seed_params(2), c.y + seed_params(1), 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
    
    xlim([seed_params(2)-3*seed_params(3) seed_params(2)+3*seed_params(3)])
    ylim([seed_params(1)-3*seed_params(3) seed_params(1)+3*seed_params(3)])
    
%     % contrast images * prf
%     subplot(3,length(im_nrs),2*length(im_nrs)+kk)
%     imEnergyPrf = imEnergyMean(im_nrs(kk),:)'.*gaufun1(seed_params(1:3));
%     imagesc(reshape(imEnergyPrf,res,res),[0 max_prf_images]); 
%     axis image
%     colorbar
end

%%% Bar plot broadband
% subplot(3,length(im_nrs),length(im_nrs)*2+[1 2]),hold on
% bar([1:length(im_nrs)],ecog_bb(im_nrs),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
% plot([1:length(im_nrs); 1:length(im_nrs)],ecog_bb_err(:,im_nrs),'k');
% xlim([0 length(im_nrs)+1])
% set(gca,'XTick',1:length(im_nrs),'XTickLabel',im_nrs)

%%% Bar plot gamma
subplot(3,length(im_nrs),length(im_nrs)*2+[1 2]),hold on
bar([1:length(im_nrs)],ecog_g(im_nrs),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
plot([1:length(im_nrs); 1:length(im_nrs)],ecog_g_err(:,im_nrs),'k');
ylim([0 max(max(ecog_g_err(:,im_nrs)))])
xlim([0 length(im_nrs)+1])
set(gca,'XTick',1:length(im_nrs),'XTickLabel',im_nrs)

%%% Bar plot model prediction
subplot(3,length(im_nrs),length(im_nrs)*2+[3 4]),hold on
bar(1:length(im_nrs),cross_NBFestimate(im_nrs),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
ylim([0 max(max(ecog_g_err(:,im_nrs)))])
xlim([0 length(im_nrs)+1])
set(gca,'XTick',1:length(im_nrs),'XTickLabel',im_nrs)

% test_params = mean(cross_NBFparams,1);
% test_params(1) = test_params(1)+4;
% test_params(3) = test_params(3)./10;
% [~,modelfit1] = helpfit_NBFsimple(stimulus,test_params,[],[]);
% plot(1:length(im_nrs),modelfit1(im_nrs)' ,'b','LineWidth',2)


% set(gcf,'PaperPositionMode','auto')
% % print('-depsc2','-r300','-painters','-tiff',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','prfSizeEffects',...
% %         ['PrfSizeEffects_sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_' modelType]))
% print('-dpsc','-r300','-painters',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','prfSizeEffects',...
%         ['PrfSizeEffects_sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_' modelType]))
% print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','prfSizeEffects',...
%         ['PrfSizeEffects_sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_' modelType]))

    
    
%% zoom into some images because the eps in teh overview file misaligns the prf and image..
max_prf_images = 0.05;

for kk = 1:length(im_nrs)
    % images + prf
    figure('Position',[0 0 500 500])
    imagesc(double(origIm.stimuli(:,:,im_nrs(kk))),[0 max(origIm.stimuli(:))]),hold on
    % Move pRF in 135x135 to original size of 800x800 
    %   half of the screen is 60 in filtered and 400 in original
    orig_x = 400 + (seed_params(1)-res/2) * 400./60;
    orig_y = 400 + (seed_params(2)-res/2) * 400./60;
    orig_sigma = seed_params(3) * 400./60;
    axis image
    colormap gray
    
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*orig_sigma);
    plot(c.x + orig_y, c.y + orig_x, 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
    [c.x, c.y] = pol2cart(c.th, 2*ones(1,numPoints)*orig_sigma);
    plot(c.x + orig_y, c.y + orig_x, 'r:') % this is just reversed because plot and imagesc are opposite, checked this with contour
    
    axis off
    
    set(gcf,'PaperPositionMode','auto')
    print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','prfSizeEffects',...
            ['PrfSizeEffects_sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_im' int2str(im_nrs(kk))]))
    print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','prfSizeEffects',...
            ['PrfSizeEffects_sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_im' int2str(im_nrs(kk))]))
end

%% Display NBF model results in 2 stimuli for several electrodes/subjects

%%%%% Pick a subject:
subject_ind = [19 19  19  19  19  19  24 24];
electrodes = [107 108 109 115 120 121 45 46];

figure('Position',[0 0 600 800])
    
for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);
    
    analysisType = 'spectra200';
    modelType = 'NBFsimple2';

    res = sqrt(size(stimulus,2)/8);

    % Load NBF model fit:
    load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveNBF',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'xys_pix','seed_params','prf_sizes',...
            'cross_NBFparams','cross_NBFestimate')

    % Load ECoG data:
    dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % Calculate ECoG gamma power percent signal change:
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = mean(100*(10.^(resamp_parms(:,:,2)-bb_base)-1),2);
    ecog_bb_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,2),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,2),.84,2))]'-bb_base)-1);
    ecog_g = mean(100*(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1),2);
    ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);
    ylims = [min(ecog_g_err(:)) max([ecog_g_err(:)])];


    res = sqrt(size(imEnergyMean,2));  % resolution of the pre-processed stimuli
    [~,xx,yy] = makegaussian2d(res,2,2,2,2);

    % Define a Gaussian centered at the prf:
    gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
        pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

    im_nrs = [10 50];
%     im_nrs = [12 69];

%     % Plot broadband responses for these two images
%     subplot(length(electrodes),6,6*ll-5),hold on
%     bar([1:length(im_nrs)],ecog_bb(im_nrs),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
%     plot([1:length(im_nrs); 1:length(im_nrs)],ecog_bb_err(:,im_nrs),'k');
%     xlim([0 length(im_nrs)+1]),ylim([0 max(max(ecog_bb_err(:,im_nrs)))])
%     set(gca,'XTick',1:length(im_nrs),'XTickLabel',im_nrs)

    % Plot gamma responses for these two images
    subplot(length(electrodes),6,6*ll-5),hold on
    bar([1:length(im_nrs)],ecog_g(im_nrs),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    plot([1:length(im_nrs); 1:length(im_nrs)],ecog_g_err(:,im_nrs),'k');
    plot(1:length(im_nrs),cross_NBFestimate(im_nrs)' ,'r','LineWidth',2)
    xlim([0 length(im_nrs)+1]),ylim([0 max(max(ecog_g_err(:,im_nrs)))])
    set(gca,'XTick',1:length(im_nrs),'XTickLabel',im_nrs)

%     % Plot gamma model predictions for the two prfs
    subplot(length(electrodes),6,6*ll-4),hold on
    bar(1:length(im_nrs),cross_NBFestimate(im_nrs),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    xlim([0 length(im_nrs)+1]),ylim([0 max(max(ecog_g_err(:,im_nrs)))])
    set(gca,'XTick',1:length(im_nrs),'XTickLabel',im_nrs)
    
    %%%%% Plot filtered (contrast) images with prf

for kk = 1:length(im_nrs)
    % images + prf
    subplot(length(electrodes),6,6*ll-4+kk)
    imagesc(double(origIm.stimuli(:,:,im_nrs(kk))),[0 max(origIm.stimuli(:))]),hold on
    % Move pRF in 135x135 to original size of 800x800 
    %   half of the screen is 60 in filtered and 400 in original
    orig_x = 400 + (seed_params(1)-res/2) * 400./60;
    orig_y = 400 + (seed_params(2)-res/2) * 400./60;
    orig_sigma = seed_params(3) * 400./60;
    axis image
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    % 1 sd
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*orig_sigma);
    plot(c.x + orig_y, c.y + orig_x, 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
    % 2 sd
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*2*orig_sigma);
    plot(c.x + orig_y, c.y + orig_x, 'r:') % this is just reversed because plot and imagesc are opposite, checked this with contour
%     xlim([orig_y-2*orig_sigma orig_y+2*orig_sigma])
%     ylim([orig_x-2*orig_sigma orig_x+2*orig_sigma])
    title(['im ' int2str(im_nrs(kk))])
    
    % contrast images + prf
    subplot(length(electrodes),6,6*ll-2+kk)
    imagesc(reshape(imEnergyMean(im_nrs(kk),:),res,res),[0 2]);
    axis image
    colormap gray
    hold on
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*seed_params(3));
    plot(c.x + seed_params(2), c.y + seed_params(1), 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
%     xlim([seed_params(2)-2*seed_params(3) seed_params(2)+2*seed_params(3)])
%     ylim([seed_params(1)-2*seed_params(3) seed_params(1)+2*seed_params(3)])
   
end

end
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','prfSizeEffects',...
        ['PrfSizeEffects_allsub_' analysisType '_' modelType]))
print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','prfSizeEffects',...
        ['PrfSizeEffects_allsub_' analysisType '_' modelType]))
