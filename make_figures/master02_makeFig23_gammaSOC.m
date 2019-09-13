% This script will make pannels for Figures 2-4 and supplementary figures
% S6-7 from: Hermes D, Petridou N, Kay K, Winawer J. 2019 An
% image-computable model for the stimulus selectivity of gamma
% oscillations. bioRxiv doi:
% https://doi.org/10.1101/583567
%
% The SOC model was adapted from:
% Kay KN, Winawer J, Rokem A, Mezer A, Wandell BA (2013) A Two-Stage
% Cascade Model of BOLD Responses in Human Visual Cortex. PLoS Comput Biol
% 9(5): e1003079. https://doi.org/10.1371/journal.pcbi.1003079
%
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

load(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

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

%%
%% Figure 2-3: example electrodes (Fig3) + average across subjects (Fig2).

% Get average broadband across all electrodes/subjects:
%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

% initialize parameters across subjects
socParams_all = zeros(length(electrodes),6);
socCOD_all = zeros(length(electrodes),2);
SOCestimate_all = zeros(length(electrodes),86);
ecog_g_all = zeros(length(electrodes),86);
ecog_g_err_all = zeros(length(electrodes),2,86);

% Get average broadband across all electrodes/subjects:
for ll = 1:length(electrodes)

    subj = subject_ind(ll);
    elec = electrodes(ll);
    
    analysisType = 'spectra200';
    modelType = 'fitSOCgammapower';

    res = sqrt(size(imEnergyMean,2));

    % load model fit
    load(fullfile(dataDir,'derivatives','gaborFilt','SOC_gamma',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seeds','cross_SOCparams','cross_SOCestimate')

    % load ecog data:
    dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj)],'ses-01','ieeg',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % Gamma power percent signal change:
    ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);
    ecog_g_all(ll,:) = ecog_g;
    ecog_g_err_all(ll,:,:) = ecog_g_err;
       
    % get mean model parameters and plot prediction
    cross_SOCparams(cross_SOCparams(:,6)<0,6) = 0; % restrictrange at 0
    cross_SOCparams(cross_SOCparams(:,6)>1,6) = 1; % restrictrange at 1
    cross_SOCparams(:,3) = abs(cross_SOCparams(:,3)); % size>0
    % write out median model parameters across 86 leave 1 outs
    socParams_all(ll,:) = median(cross_SOCparams);

    % write out coefficient of determination
    socCOD_all(ll,1) = calccod(cross_SOCestimate,ecog_g,[],0,1); % subtract mean 
    socCOD_all(ll,2) = calccod(cross_SOCestimate,ecog_g,[],0,0); % predict mean + variance

    % write out leave 1 out estimate across electrodes 
    SOCestimate_all(ll,:) = cross_SOCestimate;
end

disp(['mean SOC model performance: COD = ' int2str(mean(socCOD_all(:,2)))])



%% Plot example electrodes (Figure 4)
example_els = [3 5 7 8 9 14];

figure('Position',[0 0 470 600])
for ll = 1:length(example_els)
    elec = electrodes(example_els(ll));
    
    %%% PLOT GAMMA POWER AND SOC FIT
    subplot(8,1,ll),hold on
    bar(ecog_g_all(example_els(ll),:),1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    plot([1:86; 1:86],squeeze(ecog_g_err_all(example_els(ll),:,:)),'k');
    plot(SOCestimate_all(example_els(ll),:)','r','LineWidth',2)
    % set ylim
    ylims = [min(min(ecog_g_err_all(example_els(ll),:,:))) max(max(ecog_g_err_all(example_els(ll),:,:)))];
    ylim(ylims(1,:))
    % plot stimulus cutoffs
    stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k = 1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    set(gca,'XTick',[])
    xlim([0 87])
    ylabel(['el' int2str(elec) ' R^2=' int2str(socCOD_all(example_els(ll),2))])

end

% % save Figure 4
% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
%         ['Figure4_' modelType]))
% print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
%         ['Figure4_' modelType]))


%% Figures S6-7: all electrodes all subjects

%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

socParams_all = zeros(length(electrodes),6);
socCOD_all = zeros(length(electrodes),2);

plot_nr = 0; 
figure_nr = 1;
figure('Position',[0 0 470 600])
for ll = 1:length(electrodes)
    plot_nr = plot_nr + 1;

    subj = subject_ind(ll);
    elec = electrodes(ll);
    
    analysisType = 'spectra200';
    modelType = 'fitSOCgammapower';

    res = sqrt(size(imEnergyMean,2));
    
    % load model fit
    load(fullfile(dataDir,'derivatives','gaborFilt','SOC_gamma',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seeds','cross_SOCparams','cross_SOCestimate')

    % load ecog data:
    dataFitName = fullfile(dataDir,'derivatives','preprocessing',...
        ['sub-' int2str(subj)],'ses-01','ieeg',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % get ecog power percent signal change
    % Gamma power percent signal change:
    ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
        squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);
    ecog_g_all(ll,:) = ecog_g;
    ecog_g_err_all(ll,:,:) = ecog_g_err;
    
    %%% PLOT BROADBAND POWER AND SOC FIT
    subplot(8,1,plot_nr),hold on
    bar(ecog_g,1,'FaceColor',[.9 .9 .9],'EdgeColor',[0 0 0]);
    plot([1:86; 1:86],ecog_g_err,'k');
    plot(cross_SOCestimate' ,'r','LineWidth',2)
    % set ylim
    ylims = [min(ecog_g_err(:)) max(ecog_g_err(:))];
    ylim(ylims(1,:))
    % plot stimulus cutoffs
    stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k = 1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    set(gca,'XTick',[])
    xlim([0 87])

    % get mean model parameters and plot prediction
    cross_SOCparams(cross_SOCparams(:,6)<0,6) = 0; % restrictrange at 0
    cross_SOCparams(cross_SOCparams(:,6)>1,6) = 1; % restrictrange at 1
    cross_SOCparams(:,3) = abs(cross_SOCparams(:,3)); % size>0
    socParams_all(ll,:) = median(cross_SOCparams);

    socCOD_all(ll,1) = calccod(cross_SOCestimate,ecog_g,[],0,1);
    socCOD_all(ll,2) = calccod(cross_SOCestimate,ecog_g,[],0,0);
    
    ylabel(['el' int2str(elec) ' R^2=' int2str(socCOD_all(ll,2))])

    if mod(ll,8)==0 && ll<length(electrodes)% save figure and make a new one every 8 electrodes
        % save the figure
        set(gcf,'PaperPositionMode','auto')
        print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['FigureS6_elset' int2str(figure_nr) '_' modelType]))
        print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['FigureS6_elset' int2str(figure_nr) '_' modelType]))

        % and make a new figure
        figure_nr = figure_nr +1;
        figure('Position',[0 0 470 600])
        % reset the subplot number
        plot_nr = 0;
        
    elseif ll==length(electrodes)% save figure after last electrode
        % save the last figure
        set(gcf,'PaperPositionMode','auto')
        print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['FigureS7_elset' int2str(figure_nr) '_' modelType]))
        print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','figures',...
                ['FigureS7_elset' int2str(figure_nr) '_' modelType]))
    end
end



