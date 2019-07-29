% This script will calculate the spectral power changes for all ECoG data
% reported in: 
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. bioRxiv doi:
% https://doi.org/10.1101/583567
%
% dhermes 2019 UMC Utrecht

clear all

% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath(genpath('/Users/dora/Documents/git/ecogBasicCode/'));
addpath(genpath('/Users/dora/Documents/m-files/knkutils'));

%% epoch 

subjects = {'19','24','1001'};

for s = 1:length(subjects)
    % subject name
    subj = subjects{s};
    dataName = dir(fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_ieeg_preproc.mat']));
    stimName = dir(fullfile(dataDir,['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_events.tsv']));
    nr_runs = length(dataName);

    for data_nr = 1:nr_runs
        %%%% load preprocessed data
        load(fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
            dataName(data_nr).name));

        %%%% load stimulus information
        stim = readtable(fullfile(dataDir,['sub-' subj],'ses-01','ieeg',...
            [stimName(data_nr).name]),...
            'FileType','text','Delimiter','\t');

        %%%% notch filter data from subject 1 and 2 at 60, 120 and 180 Hz
        if ismember(subj,{'19','24'})
            data = ecog_notch(data,srate,60);
        else
            disp('no notch filter, data are pretty noise-free')
        end

        %%%% make epochs
        onset_trial = round(stim.onset*srate);%from seconds to samples
        epoch_l = 1.5; % epoch length: -0.5:1 sec
        epoch_pre = .5; % epoch pre stimulus onset time
        data_epoch = zeros(size(data,1),length(onset_trial),epoch_l*srate);
        for elec = 1:size(data,1)
            for l = 1:length(onset_trial)
                data_epoch(elec,l,:) = ...
                    data(elec,onset_trial(l)-epoch_pre*srate+1:onset_trial(l)+(epoch_l-epoch_pre)*srate);
            end
        end
        % define t - time vector for each epoch
        t = [1:epoch_l*srate]/srate - epoch_pre;    

        %%%% calculate spectra
        fft_w = window(@hann,round(srate/5)); % window width
        fft_ov = round(srate/10); % overlap
        % do not regress ERP here, because regressing out average response with a few trials can hurt 
        reg_erp = 0; % 1 to regress erp out, 0 not to, if yes, make sure to baseline correct first    

        % stimuli epochs
        fft_t = t>0.2 & t<=.5; % time segment for spectrum
        [f,data_epoch_spectra] = ...
            ecog_spectra(data_epoch(:,stim.trial_type<87,:),stim.trial_type(stim.trial_type<87),fft_w,fft_t,fft_ov,srate,reg_erp);
        stim_epoch = stim.trial_type(stim.trial_type<87);

        % ISI
        [f,data_epoch_spectra_off] = ...
            ecog_spectra(data_epoch(:,stim.trial_type==87,:),stim.trial_type(stim.trial_type==87),fft_w,fft_t,fft_ov,srate,reg_erp);
        stim_epoch_off = stim.trial_type(stim.trial_type==87);

        disp(['done run ' int2str(data_nr)])

        % SETTINGS: fft_w = window(@hann,200), fft_ov = 100, reg_erp = 0, fft_t = t>=0 & t<.5
        saveName = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
            ['sub-' subj '_ses-01_task-soc_run-' int2str(data_nr) '_spectra200_Win200.mat']);

        save(saveName,'f','data_epoch_spectra','stim_epoch','data_epoch_spectra_off','stim_epoch_off','exclude_channels','include_channels')   

    end % end run loop

end % end subject loop

%%
%% combine runs 
%%

subjects = {'19','24','1001'};

analysisType = 'spectra200_Win200';

for s = 1:length(subjects)
    % subject name
    subj = subjects{s};
    dataName = dir(fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_ieeg_preproc.mat']));
    nr_runs = length(dataName); clear dataName
   
    % initialize empty variables
    spectra = [];       % temporary variable for epoched data
    spectra_off = [];   % temporary variable for epoched data
    stims = [];         % temporary variable for stimuli list
    runs = [];          % temporary variable for stimuli list

    for data_nr = 1:nr_runs

        dataName = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
            ['sub-' subj '_ses-01_task-soc_run-' int2str(data_nr) '_' analysisType '.mat']);

        load(dataName,'f','data_epoch_spectra','stim_epoch','data_epoch_spectra_off','stim_epoch_off','exclude_channels','include_channels')   

        % concatonate with previous run
        spectra = cat(2,spectra,data_epoch_spectra);
        spectra_off = cat(2,spectra_off,data_epoch_spectra_off);
        stims = cat(1,stims,stim_epoch);
        runs = cat(1,runs,data_nr*ones(size(stim_epoch)));
    end

    saveName = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '.mat']);

    save(saveName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   
end


%
%% Load all data together and make a powerspectrum to check
%% Also print the number of trials per condition

clear all
% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

subjects = {'19','24','1001'};
s = 3;

analysisType = 'spectra200_Win200';

subj = subjects{s};

dataName = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
    ['sub-' subj '_ses-01_task-soc_allruns_' analysisType '.mat']);

load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   

disp(['this subject had ' int2str(length(find(stims==1))) ' trials per condition'])

%% plot one channel, all conditions, to check:

% channels to plot (s1: 115 | s2: 53 54 | s3: 45 46 | s4: 37 49 50 51 52 57 58 59 60):
electrodes = [50]; 

cm = jet(max(stims));

for k = 1:length(electrodes)
    figure('Color',[1 1 1],'Position',[0 0 300 300]),hold on
    elec = electrodes(k);
    
    data_fft = squeeze(spectra(elec,:,:));
    data_fft_off = squeeze(spectra_off(elec,:,:));
    
    title(['el ' int2str(elec)])
    
    % fft
    for mm = 1:max(stims)
        plot(log10(f),log10(mean(data_fft(stims==mm,:),1)),'-','Color',cm(mm,:))
    end
    plot(log10(f),log10(mean(data_fft_off,1)),'k','LineWidth',2)
    
    fill([log10(57) log10(64) log10(64) log10(57)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    fill([log10(117) log10(124) log10(124) log10(117)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    fill([log10(177) log10(184) log10(184) log10(177)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    
    xlim([log10(25) log10(200)])
    set(gca,'XTick',[log10(10) log10(25) log10(50) log10(100) log10(200)])
    set(gca,'XTickLabel',{'10','25','50','100','200'})
    
    set(gcf,'PaperPositionMode','auto')
%     print('-dpng','-r300',['./figures/spectra/' subj '_el'  int2str(electr) '_Win200'])
%     print('-depsc','-r300',['./figures/spectra/' subj '_el'  int2str(electr) '_Win200'])
%     close all
end

%% plot one channel, and select a few conditions, to check:

% channels to plot:
electrodes = [108];

% stimuli to plot:
stims_plot = [47:50 39]; % gratings contrast
cm = winter(length(stims_plot));

for k=1:length(electrodes)
    figure('Color',[1 1 1],'Position',[0 0 300 300]),hold on
    elec=electrodes(k);
    
    data_fft=squeeze(spectra(elec,:,:));
    data_fft_off=squeeze(spectra_off(elec,:,:));
    
    title(['el ' int2str(elec)])
    
    % fft
    for s = 1:length(stims_plot)
%         plot(log10(f),log10(data_fft(stims==stims_plot(s),:)),'-','Color',cm(s,:))
        plot(log10(f),log10(mean(data_fft(stims==stims_plot(s),:),1)),'-','Color',cm(s,:),'LineWidth',2)
    end
    plot(log10(f),log10(mean(data_fft_off,1)),'k','LineWidth',2)
    
%     fill([log10(57) log10(64) log10(64) log10(57)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
%     fill([log10(117) log10(124) log10(124) log10(117)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
%     fill([log10(177) log10(184) log10(184) log10(177)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    
    xlim([log10(25) log10(200)])
    set(gca,'XTick',[log10(10) log10(25) log10(50) log10(100) log10(200)])
    set(gca,'XTickLabel',{'10','25','50','100','200'})
    
    set(gcf,'PaperPositionMode','auto')
%     print('-dpng','-r300',['./figures/spectra/' subj '_el'  int2str(electr) '_Win200'])
%     print('-depsc','-r300',['./figures/spectra/' subj '_el'  int2str(electr) '_Win200'])
%     close all
end


