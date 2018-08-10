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

load(fullfile(dataDir,'soc_bids','derivatives','stimuliGaborfilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

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
s = 1; subj = subjects(s);

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
end
electrodes = [109];

res = sqrt(size(stimulus,2));  % resolution of the pre-processed stimuli

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

    % Standard SOC parameters for visual areas:
    if v_area==1
        seed_params = [xys_pix(1) xys_pix(2) xys_pix(3)*sqrt(.18) 1 .18 .93]; 
    elseif v_area==2
        seed_params = [xys_pix(1) xys_pix(2) xys_pix(3)*sqrt(.13) 1 .13 .99]; 
    elseif v_area==3
        seed_params = [xys_pix(1) xys_pix(2) xys_pix(3)*sqrt(.12) 1 .12 .99]; 
    elseif v_area==4
        seed_params = [xys_pix(1) xys_pix(2) xys_pix(3)*sqrt(.115) 1 .115 .95]; 
    end
    
    seeds = [(1+res)/2 (1+res)/2 res/4*sqrt(0.5) 1 .5 .5];
    [resultsSpace1,modelfitSpace1] = helpfit_SOC(imEnergyMean([1:38 59:68 74:78],:),[],mean(ecog_bb([1:38 59:68 74:78],:),2),seeds);
    
    seeds = [resultsSpace1.params(1:4) .5 .5];
    boot_SOCparams = zeros(size(ecog_bb,2),size(seeds,2));
    boot_SOCR = zeros(size(ecog_bb,2),1);
    boot_SOCfit = zeros(size(ecog_bb,2),size(ecog_bb,1));
    % Fit for all bootstraps:
    for kk = 1:size(ecog_bb,2) % nr of boots
        [resultsSpace,modelfitSpace] = helpfit_SOC(imEnergyMean,[],ecog_bb(:,kk),seeds);
        % put parameters together
        boot_SOCparams(kk,:) = resultsSpace.params;
        % train performance
        boot_SOCR(kk,:) = resultsSpace.trainperformance;
        % calculate modelfit
        boot_SOCfit(kk,:) = modelfitSpace;
    end
    
%     save(fullfile(dataDir,'soc_bids','derivatives','model_output','deriveSOCbb',...
%         ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_deriveSOCbblogpower']),...
%         'resultsSpace1','modelfitSpace1','xys_pix','seed_params',...
%         'boot_SOCparams','boot_SOCR','boot_SOCfit')
end

