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

%% Load preprocessed images

load(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

% compute the square root of the sum of the squares of the outputs of
% quadrature-phase filter pairs (this is the standard complex-cell energy model).
% after this step, stimulus is images x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2));
% the square root may not be typical for complex cell energy model


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
    im_deg = rad2deg(atan(17.9./50));
    electrodes = [37 43 44 45 49 50 51 52 57 58 59 60]; % S1001
end
% electrodes = [109];

res = sqrt(size(stimulus,2)/nrOrientations);  % resolution of the pre-processed stimuli

for el = 1:length(electrodes)
    elec = electrodes(el);

    % Choose an analysis type:
    analysisType = 'spectra200';
    
    % Load ECoG data:
    dataFitName = fullfile(dataDir,['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)
    
    % Gamma power percent signal change per stimulus:
    ecog_g = 100*(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1); % Actual gamma amplitude is parm3/5 width./amplitude
    
    
    %% Fit mean model

    % Initalize outputs:
    cross_MeanPred = zeros(size(ecog_g,1),1);
        
    % Estimate gain, leave one out every time for cross validation
    for kk = 1:size(ecog_g,1) % number of stimuli, leave out kk   
        % training stimuli (kk is left out)
        trainSet = setdiff([1:size(ecog_g,1)],kk);
        
        % get the mean of training stimuli
        kkEstimate = mean(mean(ecog_g(trainSet,:),2),1);
                
        % this is the estimate for the leftout stimulus
        cross_MeanPred(kk,1) = kkEstimate;
    end

    save(fullfile(dataDir,'derivatives','gaborFilt','deriveOV',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_MeanPred']),...
        'cross_MeanPred')    
end


%% Display model results for all electrodes

%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

figure('Position',[0 0 600 600])
    
for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);

    analysisType = 'spectra200';
    modelType = 'MeanPred';

    res = sqrt(size(stimulus,2)/8);

    % Load model fit:
    load(fullfile(dataDir,'derivatives','gaborFilt','deriveOV',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'cross_MeanPred')

    % Load ecog data:
    dataFitName = fullfile(dataDir,['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % Gamma power percent signal change:
    ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);

    ylims = [min(ecog_g_err(:)) max([ecog_g_err(:)])];

    subplot(15,1,ll),hold on
    
    bar(ecog_g,1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    plot([1:86; 1:86],ecog_g_err,'k');
    plot(cross_MeanPred','r','LineWidth',2)
    % Plot stimulus cutoffs:
    stim_change = [38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k = 1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    xlim([0 87]), ylim(ylims(1,:))
    set(gca,'XTick',[])
    ylabel('gamma')

end
% 
% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',fullfile(dataDir,'derivatives','gaborFilt','deriveOV',...
%         [analysisType '_allel_' modelType]))
% print('-dpng','-r300',fullfile(dataDir,'derivatives','gaborFilt','deriveOV',...
%         [analysisType '_allel_' modelType]))
    
%%
%% Plot model performance

%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

% COD output size electrodes X 2
% column 1: mean subtracted COD
% column 2: no mean subtracted COD
nbf_cod = zeros(length(electrodes),2); 

% Loop to get COD for all electrodes:
for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);

    analysisType = 'spectra200';
    modelType = 'MeanPred';

    res = sqrt(size(stimulus,2)/8);

    % Load model fit:
    load(fullfile(dataDir,'derivatives','gaborFilt','deriveOV',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'cross_MeanPred')

    % Load ecog data:
    dataFitName = fullfile(dataDir,['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % Gamma power percent signal change:
    ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    
    % Calculate leave-one-out coefficient of determination
    % mean subtracted:
    nbf_cod(ll,1) = calccod(cross_MeanPred,ecog_g,[],0,1); 
    % no mean subtracted:
    nbf_cod(ll,2) = calccod(cross_MeanPred,ecog_g,[],0,0); 
end

figure('Position',[0 0 200 300]),hold on
bar(1,mean(nbf_cod(:,1)),'w')
plot(.9+(1:length(electrodes))/50,nbf_cod(:,1),'k.')

bar(2,mean(nbf_cod(:,2)),'w')
plot(1.9+(1:length(electrodes))/50,nbf_cod(:,2),'k.')

set(gca,'XTick',[1 2],'XTickLabel',{'COD -mean','COD'})

xlim([0 3]),ylim([-10 100])

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveNBF',...
%         ['NBF_CODcross_' analysisType '_allel_' modelType]))
% print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveNBF',...
%         ['NBF_CODcross_' analysisType '_allel_' modelType]))

