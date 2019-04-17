% This script will fit the gamma percent signal change with the orientation
% variance (OV) model.
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. bioRxiv doi:
% https://doi.org/10.1101/583567
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

%% Load preprocessed images:

load(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

% compute the square root of the sum of the squares of the outputs of
% quadrature-phase filter pairs (this is the standard complex-cell energy model).
% after this step, stimulus is images x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2));
% the square root may not be typical for complex cell energy model

% we filtered the image with 8 orientations:
nrOrientations = 8;
res = sqrt(size(stimulus,2)/nrOrientations);  % resolution of the pre-processed stimuli

%% Load ECoG data and fit

%%%%% Pick a subject:
subjects = {'19','24','1001'};

for s = 1:length(subjects) 
    subj = subjects{s};

    if isequal(subj,'19') % S1
        im_deg = rad2deg(atan(17.9./50));
        electrodes = [107 108 109 115 120 121]; 
    elseif isequal(subj,'24') % S2
        im_deg = rad2deg(atan(17.9./45));
        electrodes = [45 46];
    elseif isequal(subj,'1001') % S3
        im_deg = rad2deg(atan(20./60));
        electrodes = [49 50 52 57 58 59 60]; 
    end

    for el = 1:length(electrodes)
        elec = electrodes(el);

        % get visual area: 
        [v_area] = subj_elec_visualarea(subj,elec);

        % Choose an analysis type:
        analysisType = 'spectra200';
        
        % Load ECoG data: resamp_parms for 100 bootstraps:
        dataFitName = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
            ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
        load(dataFitName)

        % Gamma power percent signal change per stimulus:
        bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
        ecog_bb = 100*(10.^(resamp_parms(:,:,2)-bb_base)-1);
        ecog_g = 100*(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1); % Actual gamma amplitude is parm3/5 width./amplitude

        %% Load SOC model results to get x and y (and sigma/sqrt(n))

        load(fullfile(dataDir,'derivatives','gaborFilt','SOC_broadband',...
            ['sub' subj '_el' int2str(elec) '_' analysisType '_fitSOCbbpower']),...
            'seeds','cross_SOCparams','cross_SOCestimate')
        
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
        [prf_s] = xy2prfsize(xys_deg,v_area);
        prf_s_pix = res.*(prf_s./im_deg);

        %% Fit gamma model

        % Exponents
        ov_exponents = [.1:.1:1];

        % Initalize outputs:
        cross_OVparams = zeros(size(ecog_bb,1),5,length(ov_exponents));
        cross_OVestimate = zeros(size(ecog_bb,1),1,length(ov_exponents));
        train_OVperformance = zeros(size(ecog_bb,1),length(ov_exponents));

        % Estimate gain, leave one out every time for cross validation
        for ii = 1:length(ov_exponents)
            % derive size from eccentricity
            seed_params = [medianParams(1:2) prf_s_pix 1 ov_exponents(ii)];    

            for kk = 1:size(ecog_bb,1) % number of stimuli, leave out kk   
                % training stimuli (kk is left out)
                trainSet = setdiff([1:size(ecog_bb,1)],kk);

                % run model through training stimuli
                [~,modelfit] = helpfit_OVexp(stimulus(trainSet,:),seed_params,[],[]);
                % get the gain
                [B] = regress(mean(ecog_g(trainSet,:),2),modelfit);

                % get training performance
                train_OVperformance(kk,ii) = calccod(B*modelfit,mean(ecog_g(trainSet,:),2),[],0,0);

                % put estimated model parameters in matrix to save
                cross_OVparams(kk,:,ii) = [seed_params(1:3) B seed_params(5)];

                % run through left out data with specific gain to get fit for
                % leftout stimulus
                [~,kkEstimate] = helpfit_OVexp(stimulus(kk,:),[seed_params(1:3) B seed_params(5)],[],[]);
                cross_OVestimate(kk,1,ii) = kkEstimate;
            end
        end

        save(fullfile(dataDir,'derivatives','gaborFilt','OV_gamma',...
            ['sub' subj '_el' int2str(elec) '_' analysisType '_OVmodel']),...
            'seed_params','cross_OVparams','cross_OVestimate','ov_exponents','train_OVperformance')
        
        disp(['Done fitting OV model for subj ' subj ' el ' int2str(elec)])

    end
end
 

%%
%% Plot model performance

%%%%% Loop over electrodes and subjects:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3

% COD output size electrodes X 2
% column 1: no mean subtracted COD
% column 2: mean subtracted COD
model_cod = zeros(length(electrodes),2,10); 
% column 1: no mean subtracted COD
train_cod = zeros(length(electrodes),1,10); 

% Loop to get COD for all electrodes:
for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);

    analysisType = 'spectra200';
    modelType = 'OVmodel';

    %%%%% Load model fit:
    load(fullfile(dataDir,'derivatives','gaborFilt','OV_gamma',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'seed_params','cross_OVparams','cross_OVestimate','ov_exponents','train_OVperformance')

    %%%%% Load ecog data:
    dataFitName = fullfile(dataDir,'derivatives','preprocessing',['sub-' int2str(subj)],'ses-01','ieeg',...
        ['sub-' int2str(subj) '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % Gamma power percent signal change:
    ecog_g = 100*(mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2));
    
    % Calculate coefficient of determination
    for ii = 1:length(ov_exponents)
        % Cross validated leave-one-out performance
        % mean subtracted:
        model_cod(ll,1,ii) = calccod(cross_OVestimate(:,:,ii),ecog_g,[],0,0); 
        % no mean subtracted:
        model_cod(ll,2,ii) = calccod(cross_OVestimate(:,:,ii),ecog_g,[],0,1); 
        
        % Calculate mean training performance 
        train_cod(ll,1,ii) = mean(train_OVperformance(:,ii));
    end
    
end

figure('Position',[0 0 300 500])
for ii = 1:length(ov_exponents)
    subplot(2,1,1),hold on
    bar(ii,mean(model_cod(:,1,ii)),'w')
    plot(ii+(1:length(electrodes))/50,model_cod(:,1,ii),'k.')
    set(gca,'XTick',1:size(cross_OVestimate,3),'XTickLabel',ov_exponents)
    ylabel('COD')
    xlim([0 size(cross_OVestimate,3)+1]),%ylim([0 100])  
    
    xlabel('n in variance^n')
end
disp(['mean OV model performance at n=.5: ' num2str(mean(model_cod(:,1,ov_exponents==.5)))])

% check training performance
figure('Position',[0 0 300 500])
for ii = 1:length(ov_exponents)
    subplot(2,1,1),hold on
    bar(ii,mean(train_cod(:,1,ii)),'w')
    plot(ii+(1:length(electrodes))/50,train_cod(:,1,ii),'k.')
    set(gca,'XTick',1:length(ov_exponents),'XTickLabel',ov_exponents)
    ylabel('train performance COD')
    xlim([0 length(ov_exponents)+1]),%ylim([0 100])  

    xlabel('n in variance^n')
end



