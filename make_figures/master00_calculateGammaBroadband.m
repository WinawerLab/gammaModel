clear all

% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath('/Users/dora/Documents/git/ecogBasicCode/')
addpath(genpath('/Users/dora/Documents/m-files/knkutils'));

%%
%% Fit gamma/bb for the electrodes in V1/2/3 with pRF within stimulus
%%

subjects = {'19','24','1001'};
electrodes_allsubjects = {[107 108 109 115 120 121], ... % S1
    [45 46], ... % S2
    [49 50 52 57 58 59 60]}; % S3
analysisType = 'spectra200';

% nr of bootstrap for resampling per stimulus
nr_boots = 100;

for s = 1:length(subjects)

    subj = subjects{s};
    electrodes = electrodes_allsubjects{s}; %

    % load data
    dataName = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '.mat']);
    load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   

    for jj = 1:length(electrodes)

        elec = electrodes(jj);
        % define output:
        resamp_parms = NaN(max(stims),nr_boots,7);

        data_fft = squeeze(spectra(elec,:,:));
        data_fft_off=squeeze(spectra_off(elec,:,:));

        f_use4fit = [30:57 65:115 126:175 186:200];
        f_sel = ismember(f,f_use4fit);
        f_alpha = find(f>=8 & f<=13);

        for k = 1:max(stims)
            disp(['fitting stimulus ' int2str(k) ' of ' int2str(max(stims))])
            % for the baseline, do not resamplel just average across 1290 trials
            data_base = mean(data_fft_off,1); % baseline

            for ii = 1:nr_boots
                % get stimulus data       
                data_fit = data_fft(stims==k,:); 

                % from baseline data, random sample with replacement
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
        end

        [a,b] = fileparts(dataName);
        dataFitName = [a '/' b '_fitEl' int2str(elec)];
        clear a b

        save(dataFitName,'resamp_parms')

    end
    clear resamp_parms spectra
end


