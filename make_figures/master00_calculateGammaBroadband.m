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



%% plot BB/G power for all stimulus conditions in a row

%%%%% Pick a subject:
subjects = {'19','24','1001'};
s = 3; subj = subjects{s};

%%%%% Pick an electrode:
electrodes_allsubjects = {[107 108 109 115 120 121], ... % S1
    [45 46], ... % S2
    [49 50 52 57 58 59 60]}; % S3
elec = 60;

% Choose an analysis type:
analysisType = 'spectra200';

% load ECoG data:
dataFitName = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
    ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
load(dataFitName,'resamp_parms')

bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
ecog_bb = 100*(10.^(squeeze(median(resamp_parms(:,:,2),2))-bb_base)-1);
ecog_bb_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,2),.16,2)) ...
    squeeze(quantile(resamp_parms(:,:,2),.84,2))]'-bb_base)-1);

ecog_g = 100*(10.^(squeeze(median(resamp_parms(:,:,3)./resamp_parms(:,:,5),2)))-1);
ecog_g_err = 100*(10.^([squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.16,2)) ...
    squeeze(quantile(resamp_parms(:,:,3)./resamp_parms(:,:,5),.84,2))]')-1);

figure('Position',[0 0 1000 100])

ylims = [...
    [min(ecog_bb_err(:)) max(ecog_bb_err(:))];...
    [0 max(ecog_g_err(:))]];

subplot(1,2,1),hold on
bar(ecog_bb,1,'b','EdgeColor',[0 0 0]);
plot([1:86; 1:86],ecog_bb_err,'k');
ylim(ylims(1,:))
ylabel('bb %change')

subplot(1,2,2),hold on
bar(ecog_g,1,'b','EdgeColor',[0 0 0]);
plot([1:86; 1:86],ecog_g_err,'k');
ylim(ylims(2,:))
ylabel('gamma %change')

for kk = 1:2
    subplot(1,2,kk),hold on
    % plot stimulus cutoffs
    stim_change = [38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for mm = 1:length(stim_change)
        plot([stim_change(mm) stim_change(mm)],ylims(kk,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    % text([19 40 46 52 55 61 70 74 79 84],zeros(10,1)+ylims(s,2)+ylims(s,2)./10,{'space','orie','grat','pl','circ','zcon','sp','zsp','coh','nm'})
    xlim([0 87])
    title(['elec ' int2str(elec)])
end

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','spectra','fit',...
%         ['BbGammaPC_' analysisType '_sub' int2str(s) '_elec' int2str(elec)]))
% print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','spectra','fit',...
%         ['BbGammaPC_' analysisType '_sub' int2str(s) '_elec' int2str(elec)]))

