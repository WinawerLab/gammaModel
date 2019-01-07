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
% the square root may not be typical for complex cell energy model

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
% 
% % sum across orientation.  after this step, stimulus is images x positions.
% imEnergyMean = blob(stimulus,2,8);


%% Load ECoG data and fit

%%%%% Pick a subject:
subjects = [19,23,24,1001];
s = 4; subj = subjects(s);

nrOrientations = 8;

if subj == 9
    im_deg = rad2deg(atan(20.7./61));
elseif subj == 19
    im_deg = rad2deg(atan(17.9./50));
    electrodes = [107 108 109 115 120 121]; % S1
elseif subj==23
    im_deg = rad2deg(atan(17.9./36));
    % electrodes = [53 54]; % S2
elseif subj==24
    im_deg = rad2deg(atan(17.9./45));
    electrodes = [45 46]; % S3
elseif subj==1001
    im_deg = rad2deg(atan(20./60));
    electrodes = [37 49 50 51 52 57 58 59 60]; % S1001 V1, V2, V3
end
% electrodes = [109];

res = sqrt(size(stimulus,2)/nrOrientations);  % resolution of the pre-processed stimuli

for el = 1:length(electrodes)
    elec = electrodes(el);
    
    % get prf from bar task, not necesary
%     % [v_area,xys,roi_labels] = subj_prf_info(subj,elec);
%     % Convert xys from degrees to pixels:
%     xys_pix = [res./2 res./2 0] + res.*(xys./im_deg);
%     xys_pix(1:2) = [res-xys_pix(2) xys_pix(1)]; % make sure that it matches images 
    % get visual area: 
    [v_area] = subj_prf_info(subj,elec);

    % Choose an analysis type:
%     analysisType = 'spectraRERP500';
    % analysisType = 'spectra';
%     analysisType = 'spectra500';
    analysisType = 'spectra200';
    
    % Load ECoG data:
    dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)
    
    % Gamma power percent signal change per stimulus:
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = 100*(10.^(resamp_parms(:,:,2)-bb_base)-1);
    ecog_g = 100*(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1); % Actual gamma amplitude is parm3/5 width./amplitude
    
    %% Load SOC model results to get x and y (and sigma/sqrt(n))
    
    modelType = 'fitSOCbbpower2';
    load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','fitSOCbb',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seeds',...
        'cross_SOCparams','cross_SOCestimate')
    % take the median of the fitting parameters, is this good?
    cross_SOCparams(cross_SOCparams(:,6)<0,6) = 0; % restrictrange c min at 0
    cross_SOCparams(cross_SOCparams(:,6)>1,6) = 1; % restrictrange c max at 1
    cross_SOCparams(:,3) = abs(cross_SOCparams(:,3)); % size>0
    medianParams = median(cross_SOCparams,1);
    
    %%%%% Derive pRF size from eccentricity %%%%%
    % Convert xys from pixels to degrees
    xys_deg = (medianParams(1:2)-[res./2 res./2])*im_deg./res;
    xys_deg(1:2) = [xys_deg(2) xys_deg(1)]; % make sure that it matches images

    % Derive pRF size from eccentricity:
    [prf_s] = xy2prfsize(xys,v_area);
    prf_s_pix = res.*(prf_s./im_deg);

    %% Fit gamma model.

    % Exponents
    ov_exponents = [.1:.1:1];

    % Initalize outputs:
    cross_OVparams = zeros(size(ecog_bb,1),5,length(ov_exponents));
    cross_OVestimate = zeros(size(ecog_bb,1),1,length(ov_exponents));
    
    
    % Estimate gain, leave one out every time for cross validation
    for ii = 1:length(ov_exponents)
        % Set the seed parameters:
        %use sigma./sqrt(n) for size
        % seed_params = [medianParams(1:2) medianParams(3)./sqrt(medianParams(5)) 1];
        % derive size from eccentricity
        seed_params = [medianParams(1:2) prf_s_pix 1 ov_exponents(ii)];    

        for kk = 1:size(ecog_bb,1) % number of stimuli, leave out kk   
            % training stimuli (kk is left out)
            trainSet = setdiff([1:size(ecog_bb,1)],kk);

            % run model through training stimuli
            [~,modelfit] = helpfit_OVexp(stimulus(trainSet,:),seed_params,[],[]);
            % get the gain
            B = regress(mean(ecog_g(trainSet,:),2),modelfit);

            % put estimated model parameters in matrix to save
            cross_OVparams(kk,:,ii) = [seed_params(1:3) B seed_params(5)];

            % run through left out data with specific gain to get fit for
            % leftout stimulus
            [~,kkEstimate] = helpfit_OVexp(stimulus(kk,:),[seed_params(1:3) B seed_params(5)],[],[]);
            cross_OVestimate(kk,1,ii) = kkEstimate;
        end
    end

    save(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_OVsimple']),...
        'seed_params',...
        'cross_OVparams','cross_OVestimate')
    % OVexp uses derived prf size from eccentricity as seed for size
    
end


%% Display model results for one electrode

%%%%% Pick a subject:
subjects = [19,23,24,1001];
s = 4; subj = subjects(s);

% %%%% Pick an electrode:
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3
% electrodes = [37 49 50 51 52 57 58 59 60]; % S1001 V1, V2, V3

analysisType = 'spectra200';
modelType = 'OVsimple';

elec = 37;
res = sqrt(size(stimulus,2)/8);

% load model fit
load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seed_params',...
        'cross_OVparams','cross_OVestimate')

% load ecog data:
dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
    'ses-01','derivatives','ieeg',...
    ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
load(dataFitName)

% ecog power
bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
ecog_bb = mean(100*(10.^(resamp_parms(:,:,2)-bb_base)-1),2);
ecog_g = mean(100*(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1),2);
ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
    squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);

ylims = [min(ecog_g_err(:)) max([ecog_g_err(:)])];

figure('Position',[0 0 1000 100])

subplot(1,2,1),hold on
bar(ecog_g,1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
plot([1:86; 1:86],ecog_g_err,'k');
% plot(cross_OVestimate' ,'r','LineWidth',2)
plot(squeeze(cross_OVestimate))
% plot stimulus cutoffs
stim_change = [38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
for k = 1:length(stim_change)
    plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
end
xlim([0 87]), ylim(ylims(1,:))
set(gca,'XTick',[1:86],'XTickLabel',[],'YTick',[0:ceil(max(ecog_g(:))/4):max(ecog_g(:))])
ylabel('gamma')

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
        ['sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_' modelType]))
print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
        ['sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_' modelType]))

%% Display model results for all electrodes

ov_exponents = [.1:.1:1];

%%%%% Pick a subject:
subject_ind = [19 19  19  19  19  19  24 24];
electrodes = [107 108 109 115 120 121 45 46];
% subject_ind = [1001 1001 1001 1001 1001 1001 1001 1001 1001];
% electrodes = [37 49 50 51 52 57 58 59 60];

figure('Position',[0 0 600 700])
    
for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);

    analysisType = 'spectra200';
    modelType = 'OVsimple';

    res = sqrt(size(stimulus,2)/8);

    % Load model fit:
    load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'seed_params',...
            'cross_OVparams','cross_OVestimate')

    % Load ecog data:
    dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % Gamma power percent signal change:
    ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);

    ylims = [min(ecog_g_err(:)) max([ecog_g_err(:)])];

    subplot(9,5,5*ll-4:5*ll-1),hold on
    
    bar(ecog_g,1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    plot([1:86; 1:86],ecog_g_err,'k');
%     plot(cross_NBFestimate','r','LineWidth',2)
    plot(squeeze(cross_OVestimate(:,:,ov_exponents==.5)),'r','LineWidth',2)
    % Plot stimulus cutoffs:
    stim_change = [38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k = 1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    xlim([0 87]), ylim(ylims(1,:))
    set(gca,'XTick',[])
    ylabel('gamma')

    %%% LOOK AT WHERE THE GAUSSIAN IS
    subplot(9,5,5*ll)
    [~,xx,yy] = makegaussian2d(res,2,2,2,2);
    imagesc(ones(size(xx)),[0 1]);
    axis image, hold on, colormap gray
    plot([res/2 res/2],[1 res],'k'),plot([1 res],[res/2 res/2],'k')
    %%% plot prf from bar/CSS model
    % gau = makegaussian2d(res,xys_pix(1),xys_pix(2),xys_pix(3),xys_pix(3),xx,yy,0,0);
    % imagesc(gau);
    % axis image, hold on, colorbar
    % look at the prf from the SOC fit:
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    % plot 1 sd
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*seed_params(3));
    plot(c.x + seed_params(2), c.y + seed_params(1), 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
    % plot 2 sd
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*2*seed_params(3));
    plot(c.x + seed_params(2), c.y + seed_params(1), 'r:') % this is just reversed because plot and imagesc are opposite, checked this with contour
    axis off
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
        [analysisType '_allel_' modelType]))
print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
        [analysisType '_allel_' modelType]))
% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
%         [analysisType '_allel_' modelType '_SubChaam']))
% print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
%         [analysisType '_allel_' modelType '_SubChaam']))
    
%%
%% Plot model performance
% Exponents
ov_exponents = [.1:.1:1];

%%%%% Loop over electrodes and subjects:
subject_ind = [19 19  19  19  19  19  24 24];
electrodes = [107 108 109 115 120 121 45 46];
% subject_ind = [1001 1001 1001 1001 1001 1001 1001 1001 1001];
% electrodes = [37 49 50 51 52 57 58 59 60];
subject_ind = [19 19  19  19  19  19  24 24 1001 1001 1001 1001 1001 1001 1001 1001 1001];
electrodes = [107 108 109 115 120 121 45 46 37 49 50 51 52 57 58 59 60];

% COD output size electrodes X 2
% column 1: mean subtracted COD
% column 2: no mean subtracted COD
nbf_cod = zeros(length(electrodes),2); 

% Loop to get COD for all electrodes:
for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);

    analysisType = 'spectra200';
    modelType = 'OVsimple';

    res = sqrt(size(stimulus,2)/8);

    % Load model fit:
    load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'seed_params',...
            'cross_OVparams','cross_OVestimate')

    % Load ecog data:
    dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % Gamma power percent signal change:
    ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    
    % Calculate leave-one-out coefficient of determination
    for ii = 1:size(cross_OVestimate,3)
        % mean subtracted:
        nbf_cod(ll,1,ii) = calccod(cross_OVestimate(:,:,ii),ecog_g,[],0,1); 
        % no mean subtracted:
        nbf_cod(ll,2,ii) = calccod(cross_OVestimate(:,:,ii),ecog_g,[],0,0); 
    end
end
%%
figure('Position',[0 0 300 500])
for ii = 1:size(cross_OVestimate,3)
    subplot(2,1,1),hold on
    bar(ii,mean(nbf_cod(:,1,ii)),'w')
    plot(ii+(1:length(electrodes))/50,nbf_cod(:,1,ii),'k.')
    set(gca,'XTick',1:size(cross_OVestimate,3),'XTickLabel',ov_exponents)
    ylabel('COD-mean')
    xlim([0 size(cross_OVestimate,3)+1]),ylim([0 100])  
    
    subplot(2,1,2),hold on
    bar(ii,mean(nbf_cod(:,2,ii)),'w')
    plot(ii+(1:length(electrodes))/50,nbf_cod(:,2,ii),'k.')
    set(gca,'XTick',1:size(cross_OVestimate,3),'XTickLabel',ov_exponents)
    ylabel('COD')
    xlim([0 size(cross_OVestimate,3)+1]),ylim([0 100])
    xlabel('n in variance^n')
end
% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
%         ['OV_CODcross_' analysisType '_allel_' modelType]))
% print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveOV',...
%         ['OV_CODcross_' analysisType '_allel_' modelType]))

