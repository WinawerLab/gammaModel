clear

% set paths:
gammaModelCodePath;
% put here the path to the data, and make sure that line 29 loads the data
% correctly:
dataDir = gammaModelDataPath;

%%
%% Fit gamma/bb for the electrodes in V1/2/3 with pRF within stimulus
%% STIMULUS CONDITIONS
%%

% subjects = {'19','24','1001'};
% electrodes_allsubjects = {[107 108 109 115 120 121], ... % S1
%    [45 46], ... % S2
%    [49 50 52 57 58 59 60]}; % S3

subjects                = [ 19  19  19  19  19  19 24 24 1001 1001 1001 1001 1001 1001 1001];
electrodes_allsubjects  = [107 108 109 115 120 121 45 45   49   50   52   57   58   59   60];

analysisType = 'spectra200';


hpc_job_number = str2num(getenv('SLURM_ARRAY_TASK_ID'));


% nr of bootstrap for resampling per stimulus
nr_boots = 1000;

subj = num2str(subjects(hpc_job_number));
elec = electrodes_allsubjects(hpc_job_number);

% load data
dataName = fullfile(dataDir,...
    ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '.mat']);
load(dataName,'f','spectra','spectra_off','stims')


[a,b] = fileparts(dataName);
dataFitName = [a '/' b '_fitEl' int2str(elec)];

% define output:
resamp_parms = NaN(max(stims),nr_boots,7);

data_fft = squeeze(spectra(elec,:,:));
data_fft_off = squeeze(spectra_off(elec,:,:));

f_use4fit = [30:57 65:115 126:175 186:200];
f_sel = ismember(f,f_use4fit);
f_alpha = find(f>=8 & f<=13);

for k = 1:max(stims)
    disp(['fitting stimulus ' int2str(k) ' of ' int2str(max(stims))])
    % for the baseline, do not resample just average across the trials
    data_base = mean(data_fft_off,1); % baseline
    
    for ii = 1:nr_boots
        % get stimulus data
        data_fit = data_fft(stims==k,:);
        
        % from stimulus data, random sample with replacement
        trial_set = randsample(size(data_fit,1),size(data_fit,1),true);
        
        % average across resampled trials
        data_fit = mean(data_fit(trial_set,:),1);
        
        % do the fitting
        [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
            ecog_fitgamma(f,f_use4fit,data_base,data_fit);
        resamp_parms(k,ii,1) = out_exp(1); % this is the slope used in all cases
        resamp_parms(k,ii,2) = bb_amp;
        resamp_parms(k,ii,3) = gamma_amp;
        resamp_parms(k,ii,4) = gamma_freq;
        resamp_parms(k,ii,5) = gamma_width;
        resamp_parms(k,ii,6) = out_exp(2); % this is the baseline intercept
        clear trial_set
        
        % calculate alpha change
        resamp_parms(k,ii,7) = mean(log10(data_fit(f_alpha)) - log10(data_base(f_alpha)));
    end
    
    save(dataFitName,'resamp_parms')

end


clear resamp_parms spectra



%%
%% Fit gamma/bb for the electrodes in V1/2/3 with pRF within stimulus
%% BASELINE
%%

% subjects = {'19','24','1001'};
% electrodes_allsubjects = {[107 108 109 115 120 121], ... % S1
%     [45 46], ... % S2
%     [49 50 52 57 58 59 60]}; % S3

analysisType = 'spectra200';


% load data
dataName = fullfile(dataDir,...
    ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '.mat']);
load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')

[a,b] = fileparts(dataName);
dataFitName = [a '/' b '_fitBaseEl' int2str(elec)];
clear a b


% define output:
resamp_parms_baseline = NaN(1,nr_boots,7);

data_fft = squeeze(spectra(elec,:,:));
data_fft_off = squeeze(spectra_off(elec,:,:));

f_use4fit = [30:57 65:115 126:175 186:200];
f_sel = ismember(f,f_use4fit);
f_alpha = find(f>=8 & f<=13);

% for the baseline, do not resample just average across the trials
data_base = mean(data_fft_off,1); % baseline

for ii = 1:nr_boots
    % get baseline data to estimate
    data_fit = data_fft_off;
    
    % from stimulus data, random sample with replacement
    n_stims_draw = length(find(stims==1));
    trial_set = randsample(size(data_fit,1),n_stims_draw,true);
    
    % average across resampled trials
    data_fit = mean(data_fit(trial_set,:),1);
    
    % do the fitting
    [out_exp,bb_amp,gamma_amp,gamma_freq,gamma_width,fit_f2] = ...
        ecog_fitgamma(f,f_use4fit,data_base,data_fit);
    resamp_parms_baseline(1,ii,1) = out_exp(1); % this is the slope used in all cases
    resamp_parms_baseline(1,ii,2) = bb_amp;
    resamp_parms_baseline(1,ii,3) = gamma_amp;
    resamp_parms_baseline(1,ii,4) = gamma_freq;
    resamp_parms_baseline(1,ii,5) = gamma_width;
    resamp_parms_baseline(1,ii,6) = out_exp(2); % this is the baseline intercept
    clear trial_set
    
    % calculate alpha change
    resamp_parms_baseline(1,ii,7) = mean(log10(data_fit(f_alpha)) - log10(data_base(f_alpha)));
    
end

save(dataFitName,'resamp_parms_baseline')


clear resamp_parms_baseline spectra

