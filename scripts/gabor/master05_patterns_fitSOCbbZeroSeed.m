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

%% Load preprocessed images and divisive normalization:

load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

% compute the square root of the sum of the squares of the outputs of
% quadrature-phase filter pairs (this is the standard complex-cell energy model).
% after this step, stimulus is images x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2));

% compute the population term in the divisive-normalization equation.
% this term is simply the average across the complex-cell outputs
% at each position (averaging across orientation).
stimulusPOP = blob(stimulus,2,8)/8;

% repeat the population term for each of the orientations
stimulusPOP = upsamplematrix(stimulusPOP,8,2,[],'nearest');

% apply divisive normalization to the complex-cell outputs.  there are two parameters
% that influence this operation: an exponent term (r) and a semi-saturation term (s).
% the parameter values specified here were determined through a separate fitting
% procedure (see paper for details).  for the purposes of this script, we will
% simply hard-code the parameter values here and not worry about attempting to fit
% the parameters.
r = 1;
s = 0.5;
stimulus = stimulus.^r ./ (s.^r + stimulusPOP.^r);
clear stimulusPOP;

% sum across orientation.  after this step, stimulus is images x positions.
imEnergyMean = blob(stimulus,2,8);

%% Load ECoG data and fit

%%%%% Pick a subject:
subjects = [19,23,24];
for s = [1 3]; 
subj = subjects(s);

if subj == 9
    im_deg = rad2deg(atan(20.7./61));
elseif subj == 19
    im_deg = rad2deg(atan(17.9./50));
%     electrodes = [107 108 109 115 120 121]; % S1
    electrodes = [108 109 115 120 121]; % S1
elseif subj==23
    im_deg = rad2deg(atan(17.9./36));
    % electrodes = [53 54]; % S2
elseif subj==24
    im_deg = rad2deg(atan(17.9./45));
%     electrodes = [45 46]; % S3
    electrodes = [46];
end
% electrodes = [45];

res = sqrt(size(imEnergyMean,2));  % resolution of the pre-processed stimuli

for el = 1:length(electrodes)
    elec = electrodes(el);

    [v_area,xys,roi_labels] = subj_prf_info(subj,elec);
    % Convert xys from degrees to pixels
    xys_pix = [res./2 res./2 0] + res.*(xys./im_deg);
    xys_pix(1:2) = [res-xys_pix(2) xys_pix(1)]; % make sure that it matches images
    
    % Choose an analysis type:
%     analysisType = 'spectraRERP500';
    % analysisType = 'spectra';
%     analysisType = 'spectra500';
    analysisType = 'spectra200';

%     % load resamp_parms_even and resamp_parms_odd for even/odd:
%     dataFitName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '_evenodd.mat'];
%     load(dataFitName)
    % load resamp_parms for 100 bootstraps:
    dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)
    
    % Broadband estimate, one value per image:
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = 10.^(resamp_parms(:,:,2)-bb_base)-1;
%     ecog_bb = resamp_parms(:,:,2)-bb_base;
    % ecog_g = 10.^resamp_parms(:,:,3)-1;
    
    %% Fit SOC model on bb to confirm prf location?

    % The SOC model parameters that we will be fitting are [R C S G N] where:
    %   R is the row index of the center of the 2D Gaussian
    %   C is the column index of the center of the 2D Gaussian
    %   S is the size of the 2D Gaussian
    %   G is a gain parameter
    %   N is an exponent

    % Get a good gain seed:
    gain_seed = max(mean(ecog_bb,2));
        
    % Seeds for fitting n and c
    Ns = [.1 .3 .5 .7 .9 1];
    Cs = [.1 .4 .7 .8 .9 .95 1];
    seeds = [];
    for p = 1:length(Ns)
      for q = 1:length(Cs)
        seeds = cat(1,seeds,[(1+res)/2 (1+res)/2 res/4*sqrt(Ns(p)) gain_seed Ns(p) Cs(q)]);  
      end
    end
        
    % Initalize outputs:
    cross_SOCparams = zeros(size(ecog_bb,1),6);
    cross_SOCestimate = zeros(size(ecog_bb,1),1);
    
    % Estimate gain, n, c, leave one out every time for cross validation
    for kk = 1:size(ecog_bb,1) % number of stimuli, leave out kk   
        
        % training stimuli (kk is left out)
        trainSet = setdiff([1:size(ecog_bb,1)],kk);
        
        % calculate decent model for space stimuli  
        [resultsSpace1] = helpfit_SOC2(...
            imEnergyMean(trainSet(ismember(trainSet,1:38)),:),[],...
            mean(ecog_bb(trainSet(ismember(trainSet,1:38)),:),2),seeds);
        % update seeds with the decent model
        seeds(:,1:3) = repmat(resultsSpace1.params(1:3),size(seeds,1),1);
        
        % calculate model
        [resultsSpace,modelfitSpace] = helpfit_SOC2(...
            imEnergyMean(trainSet,:),[],...
            mean(ecog_bb(trainSet,:),2),seeds);
        
        % train parameters
        cross_SOCparams(kk,:) = resultsSpace.params;

        % estimate for the left-out stimulus kk
        [~,kkEstimate] = helpfit_SOC2(imEnergyMean(kk,:),resultsSpace.params,[],[]);
        cross_SOCestimate(kk) = kkEstimate;
        clear kkEstimate
        
    end

    save(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','fitSOCbb',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_fitSOCbbpower2']),...
        'xys_pix','seeds',...
        'cross_SOCparams','cross_SOCestimate')
end
end



%% 
%% Display results for one electrode

%%%%% Pick a subject:
subjects = [19,23,24];
s = 3; subj = subjects(s);

% %%%% Pick an electrode:
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3

analysisType = 'spectra200';
modelType = 'fitSOCbbpower2';

elec = 46;
res = 240;

% load model fit
load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','fitSOCbb',...
    ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
    'xys_pix','seeds',...
    'cross_SOCparams','cross_SOCestimate');

% load ecog data:
dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
    'ses-01','derivatives','ieeg',...
    ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
load(dataFitName)

% ecog power
bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
ecog_bb = mean(10.^(resamp_parms(:,:,2)-bb_base)-1,2);
ecog_g = 10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1;

figure('Position',[0 0 1200 160])

ylims = [min(ecog_bb(:)) max([ecog_bb(:); cross_SOCestimate(:)+.1])];

subplot(1,3,1:2),hold on
bar(ecog_bb,1,'b','EdgeColor',[0 0 0]);
plot(cross_SOCestimate' ,'g','LineWidth',2)
title(['elec ' int2str(elec)])
% plot stimulus cutoffs
stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
for k=1:length(stim_change)
    plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
end
xlim([0 87]), ylim(ylims(1,:))
set(gca,'YTick',[0:1:floor(max(ecog_bb))])
ylabel('bb')

subplot(1,3,3)
%%% LOOK AT WHERE THE GAUSSIAN IS
[~,xx,yy] = makegaussian2d(res,2,2,2,2);
% gau = makegaussian2d(res,results.params(1),results.params(2),results.params(3)/sqrt(results.params(5)),results.params(3)/sqrt(results.params(5)),xx,yy,0,0);
gau = makegaussian2d(res,xys_pix(1),xys_pix(2),xys_pix(3),xys_pix(3),xx,yy,0,0);
imagesc(gau);
colorbar
axis image
hold on
% look at the prf from the SOC fit:
numPoints = 50;
c.th = linspace(0,2*pi, numPoints);
for kk = 1:size(cross_SOCparams,1)
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*cross_SOCparams(kk,3)./sqrt(cross_SOCparams(kk,5)));
    plot(c.x + cross_SOCparams(kk,2), c.y + cross_SOCparams(kk,1), 'k') % this is just reversed because plot and imagesc are opposite, checked this with contour
end
title('bar pRF (color) and SOC pRF (black)')
 
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','fitSOCbb',...
        ['sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_' modelType]))
print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','fitSOCbb',...
        ['sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_' modelType]))

%% Figure all subjects

%%%%% Pick a subject:
subject_ind = [19 19  19  19  19  19  24 24];
electrodes = [107 108 109 115 120 121 45 46];

figure('Position',[0 0 1400 800])
    
for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);
    
    analysisType = 'spectra200';
    modelType = 'fitSOCbbpower2';

    res = 240;

    % load model fit
    load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','fitSOCbb',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'xys_pix','seeds',...
        'cross_SOCparams','cross_SOCestimate')

    % load ecog data:
    dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % ecog power
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = mean(10.^(resamp_parms(:,:,2)-bb_base)-1,2);
    ecog_bb_yneg = median(10.^(resamp_parms(:,:,2)-bb_base)-1,2)-...
        quantile(10.^(resamp_parms(:,:,2)-bb_base)-1,.16,2);
    ecog_bb_ypos = quantile(10.^(resamp_parms(:,:,2)-bb_base)-1,.84,2)-...
        median(10.^(resamp_parms(:,:,2)-bb_base)-1,2);
    ecog_g = 10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1;

    ylims = [min(ecog_bb-ecog_bb_yneg) max([ecog_bb+ecog_bb_ypos; cross_SOCestimate(:)+.1])];

    subplot(8,2,2*ll-1),hold on
    bar(ecog_bb,1,'b','EdgeColor',[0 0 0]);
%     errorbar([1:length(ecog_bb)],ecog_bb,ecog_bb_yneg,ecog_bb_ypos,'k')
    plot(cross_SOCestimate' ,'g','LineWidth',2)
    ylim(ylims(1,:))
    
    % plot stimulus cutoffs
    stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k = 1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    xlim([0 87])
    set(gca,'XTick',stim_change,'XTickLabel',[],...
        'YTick',[0:2:floor(max(ecog_bb))])
    ylabel(['bb el ' int2str(elec)])

    % get mean model parameters and plot prediction
    cross_SOCparams(cross_SOCparams(:,6)<0,6) = 0; % restrictrange at 0
    cross_SOCparams(cross_SOCparams(:,6)>1,6) = 1; % restrictrange at 1
    cross_SOCparams(:,3) = abs(cross_SOCparams(:,3)); % size>0
    medianParams = median(cross_SOCparams);
    % plot median
%     [~,bbEstimate] = helpfit_SOC2(imEnergyMean,medianParams,[],[]);
%     plot(bbEstimate,'r','LineWidth',2)

    subplot(4,4,3),hold on
    bar(ll,median(cross_SOCparams(:,5)),'w')
%     plot(ll,cross_SOCparams(:,5),'k.')
    % 68% ci  ~ 1sd
    errorbar(ll,median(cross_SOCparams(:,5)),...
        median(cross_SOCparams(:,5))-quantile(cross_SOCparams(:,5),.16),...
        quantile(cross_SOCparams(:,5),.84)-median(cross_SOCparams(:,5)),'k')
    
    subplot(4,4,4),hold on
    bar(ll,median(cross_SOCparams(:,6)),'w')
%     plot(ll,cross_SOCparams(:,6),'k.')
    % 68% ci  ~ 1sd
    errorbar(ll,median(cross_SOCparams(:,6)),...
        median(cross_SOCparams(:,6))-quantile(cross_SOCparams(:,6),.16),...
        quantile(cross_SOCparams(:,6),.84)-median(cross_SOCparams(:,6)),'k')

    subplot(4,4,7),hold on
    bar(ll,calccod(cross_SOCestimate,ecog_bb,[],0,1),'w')

    subplot(4,4,8),hold on
    bar(ll,calccod(cross_SOCestimate,ecog_bb,[],0,0),'w')
end

subplot(4,4,3),hold on
xlim([0 9]),ylim([0 1.1])
title('n parameter')
xlabel('electrode')
set(gca,'XTick',[1:8],'XTickLabel',electrodes)

subplot(4,4,4),hold on
xlim([0 9]),ylim([0 1.1])
title('c parameter')
xlabel('electrode')
set(gca,'XTick',[1:8],'XTickLabel',electrodes)

subplot(4,4,7),hold on
title('cod mean subtracted')
xlim([0 9]),ylim([0 110])
set(gca,'XTick',[1:8],'XTickLabel',electrodes)

subplot(4,4,8),hold on
title('cod no mean subtracted')
xlim([0 9]),ylim([0 110])
set(gca,'XTick',[1:8],'XTickLabel',electrodes)

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','fitSOCbb',...
        [analysisType '_allel_' modelType]))
print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','fitSOCbb',...
        [analysisType '_allel_' modelType]))


