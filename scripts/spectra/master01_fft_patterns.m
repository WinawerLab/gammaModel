clear all

addpath('~/Documents/git/ecogBasicCode/')
addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);
addpath(genpath('~/Documents/m-files/knkutils'));

%% epoch 
% NOTE: later add exclusion bad epochs

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = {19,24,1001};

for s = 1:length(subjects)
% subject name
subj = subjects{s};
dataName = dir([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/sub-' int2str(subj) '_ses-01_task-soc_run-*_ieeg_preproc.mat']);
stimName = dir([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/sub-' int2str(subj) '_ses-01_task-soc_run-*_events.tsv']);
nr_runs = length(dataName);

for data_nr = 1:nr_runs
    %%%% load preprocessed data
    load([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/' dataName(data_nr).name]);
    
    %%%% load stimulus information
    stim = readtable([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/' stimName(data_nr).name],...
        'FileType','text','Delimiter','\t');

    %%%% notch filter data at 60, 120 and 180 Hz
    data = ecog_notch(data,srate);
    
    %%%% make epochs
    onset_trial = round(stim.onset*srate);%from seconds to samples
    epoch_l = 1.2; % epoch length: -0.2:1 sec
    data_epoch = zeros(size(data,2),length(onset_trial),epoch_l*srate);
    for elec = 1:size(data,2)
        for l = 1:length(onset_trial)
            data_epoch(elec,l,:) = ...
                data(onset_trial(l)-.2*srate+1:onset_trial(l)+(epoch_l-.2)*srate,elec)';
        end
    end
    % define t - time vector for each epoch
    t = [1:epoch_l*srate]/srate - 0.2;    

    %%%% baseline correct
    % data_epoch = ecog_baselinesubtract(data_epoch,t>-.1 & t<0);

    %%%% calculate spectra
    fft_w = window(@hann,200); % window width
    fft_ov = 100; % overlap
    % do not regress ERP here, because regressing out average response with a few trials can hurt 
    reg_erp = 0; % 1 to regress erp out, 0 not to, if yes, make sure to baseline correct first    
    
    % stimuli epochs
    fft_t = t>0 & t<=.5; % time segment for spectrum
    [f,data_epoch_spectra] = ...
        ecog_spectra(data_epoch(:,stim.trial_type<87,:),stim.trial_type(stim.trial_type<87),fft_w,fft_t,fft_ov,srate,reg_erp);
    stim_epoch = stim.trial_type(stim.trial_type<87);
    
    % ISI
    [f,data_epoch_spectra_off] = ...
        ecog_spectra(data_epoch(:,stim.trial_type==87,:),stim.trial_type(stim.trial_type==87),fft_w,fft_t,fft_ov,srate,reg_erp);
    stim_epoch_off = stim.trial_type(stim.trial_type==87);

    disp(['done run ' int2str(data_nr)])
    
    % Do not regress ERP out for each stimulus, better to take 200-500 ms
    % window.
    
    % SETTINGS: fft_w = window(@hann,300), fft_ov = 100, reg_erp = 0, fft_t = t>=0.2 & t<.5
%     saveName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_run-' int2str(data_nr) '_spectra.mat'];

    % SETTINGS: fft_w = window(@hann,300), fft_ov = 100, reg_erp = 1, fft_t = t>=0.2 & t<.5
%     saveName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_run-' int2str(data_nr) '_spectraRERP.mat'];

    % SETTINGS: fft_w = window(@hann,500), fft_ov = 100, reg_erp = 0, fft_t = t>=0 & t<.5
%     saveName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_run-' int2str(data_nr) '_spectra500.mat'];

    % SETTINGS: fft_w = window(@hann,500), fft_ov = 100, reg_erp = 1, fft_t = t>=0 & t<.5
%     saveName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_run-' int2str(data_nr) '_spectraRERP500.mat'];

    % SETTINGS: fft_w = window(@hann,200), fft_ov = 100, reg_erp = 0, fft_t = t>=0 & t<.5
    saveName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_run-' int2str(data_nr) '_spectra200.mat'];
    
    % SETTINGS: fft_w = window(@hann,200), fft_ov = 100, reg_erp = 1, fft_t = t>=0 & t<.5
%     saveName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_run-' int2str(data_nr) '_spectraRERP200.mat'];

    % SETTINGS: fft_w = window(@hann,100), fft_ov = 50, reg_erp = 0, fft_t = t>=0 & t<.5
%     saveName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_run-' int2str(data_nr) '_spectra100.mat'];

    % SETTINGS: fft_w = window(@hann,300), fft_ov = 150, reg_erp = 0, fft_t = t>=0 & t<.5
%     saveName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_run-' int2str(data_nr) '_spectra300.mat'];

    save(saveName,'f','data_epoch_spectra','stim_epoch','data_epoch_spectra_off','stim_epoch_off','exclude_channels','include_channels')   
    
end % end run loop

end % end subject loop

%%
%% combine runs 
%%

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [19,23,24];

% analysisType = 'spectra';
% analysisType = 'spectraRERP';
% analysisType = 'spectra500';
% analysisType = 'spectraRERP500';
% analysisType = 'spectra200';
analysisType = 'spectraRERP200';
% analysisType = 'spectra100';
analysisType = 'spectra300';

for s = 1:length(subjects)
    % subject name
    subj = subjects(s);
    dataName = dir([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/sub-' int2str(subj) '_ses-01_task-soc_run-*_ieeg_preproc.mat']);
    nr_runs = length(dataName); clear dataName
   
    % initialize empty variables
    spectra = [];  % temporary variable for epoched data
    spectra_off = [];  % temporary variable for epoched data
    stims = [];       % temporary variable for stimuli list
    runs = [];       % temporary variable for stimuli list

    for data_nr = 1:nr_runs

        dataName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_run-' int2str(data_nr) '_' analysisType '.mat'];

        load(dataName,'f','data_epoch_spectra','stim_epoch','data_epoch_spectra_off','stim_epoch_off','exclude_channels','include_channels')   

        % concatonate with previous run
        spectra = cat(2,spectra,data_epoch_spectra);
        spectra_off = cat(2,spectra_off,data_epoch_spectra_off);
        stims = cat(1,stims,stim_epoch);
        runs = cat(1,runs,data_nr*ones(size(stim_epoch)));
    end

    saveName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '.mat'];

    save(saveName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   

end

%% load all data together
%%
clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids_sourcedir';
subjects = [19,23,24,1001];
s = 3;

% analysisType = 'spectra';
% analysisType = 'spectra500';
% analysisType = 'spectra300';
analysisType = 'spectra200';
% analysisType = 'spectra100';
% analysisType = 'spectraRERP';
% analysisType = 'spectraRERP500';
% analysisType = 'spectraRERP200';

subj = subjects(s);
dataName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '.mat'];
load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   


%% plot one channel, all conditions, to check:

% channels to plot (s1: 115 | s2: 53 54 | s3: 45 46 | s4: 37 49 50 51 52 57 58 59 60):
electrodes=[58]; 

cm = jet(max(stims));

for k=1:length(electrodes)
    figure('Color',[1 1 1],'Position',[0 0 300 300]),hold on
    elec=electrodes(k);
    
    data_fft=squeeze(spectra(elec,:,:));
    data_fft_off=squeeze(spectra_off(elec,:,:));
    
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
%     print('-dpng','-r300',['./figures/spectra/' subj '_el'  int2str(electr)])
%     print('-depsc','-r300',['./figures/spectra/' subj '_el'  int2str(electr)])
%     close all
end

%% plot one channel, and select a few conditions, to check:

% channels to plot:
electrodes = [58];

% stimuli to plot:
% stims_plot = [39:46]; % gratings orientation
% stims_plot = [47:50 39]; % gratings contrast
% stims_plot = [11 29]; 
% stims_plot = [79:82]; 
stims_plot = [74:78]; 
% cm = jet(length(stims_plot));
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
%     print('-dpng','-r300',['./figures/spectra/' subj '_el'  int2str(electr)])
%     print('-depsc','-r300',['./figures/spectra/' subj '_el'  int2str(electr)])
%     close all
end


%% compare analyses: load all data together
%%
clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subjects = [19,23,24];
s = 3;
% Pick an electrode:
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3
elec = 45;

analysisTypes = {'spectra500','spectraRERP500','spectra','spectraRERP','spectra200','spectraRERP200','spectra100'};
plot_colors = [[0 0 0];[.5 .5 .5];[1 0 0];[1 .5 0];[0 0 1];[0 1 1];[1 0 1]];
subj = subjects(s);

% stims_plot = [39:46]; % gratings orientation
% stims_plot = [11 29]; % full field pattern
stims_plot = {[11 29],39:46}; 
stims_names = {'full field pattern','grating orientations'};

figure('Position',[0 0 800 300])

for aa = 1:length(analysisTypes)
    thisAnalysisType = analysisTypes{aa};

    dataName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' thisAnalysisType '.mat'];
    load(dataName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   
    
    for ss = 1:length(stims_plot)
        subplot(1,length(stims_plot)+1,ss),hold on
    
        % plot stimulus spectra for all stimuli:
        data_fft=squeeze(mean(log10(spectra(elec,ismember(stims,stims_plot{ss}),:)),2));
        plot(log10(f),data_fft,'-','Color',plot_colors(aa,:),'LineWidth',2)
    
        % plot rest spectra for all stimuli:
        data_fft_off=squeeze(mean(log10(spectra_off(elec,:,:)),2));
        plot(log10(f),data_fft_off,':','Color',plot_colors(aa,:),'LineWidth',2)
    
        % plot details in all plots:    
        fill([log10(57) log10(64) log10(64) log10(57)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
        fill([log10(117) log10(124) log10(124) log10(117)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
        fill([log10(177) log10(184) log10(184) log10(177)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
        xlim([log10(5) log10(200)])
        set(gca,'XTick',[log10(10) log10(25) log10(50) log10(100) log10(200)])
        set(gca,'XTickLabel',{'10','25','50','100','200'})
        title(['el ' int2str(elec) ' stims ' stims_names{ss}])
    end
end

subplot(1,length(stims_plot)+1,length(stims_plot)+1),hold on
for aa = 1:length(analysisTypes)
    plot(0,aa,'.','Color',plot_colors(aa,:))
    text(0.2,aa,analysisTypes{aa},'Color',plot_colors(aa,:))
end
ylim([0 length(analysisTypes)+1])

set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/analysisCompare/analysisCompare_sub-' int2str(subj) '_el'  int2str(elec)])
% print('-depsc','-r300',['./figures/analysisCompare/analysisCompare_sub-' int2str(subj) '_el'  int2str(elec)])



