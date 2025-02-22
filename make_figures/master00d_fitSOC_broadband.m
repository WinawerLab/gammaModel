% This script will fit the broadband percent signal change with the SOC
% model.
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. bioRxiv doi:
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

%% Load ECoG data and fit
% This takes a while to run: about 1-2 hours per electrode, because of the
% resampling.

%%%%% Pick a subject:
subjects = {'19','24','1001'};

for s = 1:3
    subj = subjects{s};

    if isequal(subj,'19') % S1
        im_deg = rad2deg(atan(17.9./50));
        electrodes = [107 108 109 115 120 121]; 
    elseif isequal(subj,'24') % S2
        im_deg = rad2deg(atan(17.9./45));
        electrodes = [45 46]; 
    elseif isequal(subj,'1001') % S3
        im_deg = rad2deg(atan(17.9./50));
        electrodes = [49 50 52 57 58 59 60]; 
    end

    res = sqrt(size(imEnergyMean,2));  % resolution of the pre-processed stimuli

    for el = 1:length(electrodes)
        elec = electrodes(el);

        % Choose an analysis type:
        analysisType = 'spectra200';

        % Load ECoG data: resamp_parms for 100 bootstraps
        dataFitName = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
            ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
        load(dataFitName)

        % Broadband power percent signal change, one value per image:
        bb_base = resamp_parms(1,1,6); % The baseline is the same for resamp_parms(:,:,6)
        ecog_bb = 100*(10.^(resamp_parms(:,:,2)-bb_base)-1);

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

        save(fullfile(dataDir,'derivatives','gaborFilt','SOC_broadband',...
            ['sub' subj '_el' int2str(elec) '_' analysisType '_fitSOCbbpower']),...
            'seeds','cross_SOCparams','cross_SOCestimate')
    end
end
 

