clear all

addpath('~/Documents/git/ecogBasicCode/')
addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);
addpath(genpath('/Users/dorahermes-miller/Documents/m-files/knkutils'));

%% epoch and fft

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [9];

for s = 1:length(subjects)
    % subject name
    if subjects(s)<10
        subj = ['0' int2str(subjects(s))];
    else
        subj = int2str(subjects(s));
    end

    dataName = dir([dataRootPath '/sub-' subj '/ses-01/ieeg/sub-' subj '_task-objects_run-*_ieeg_preproc.mat']);
    stimName = dir([dataRootPath '/sub-' subj '/ses-01/ieeg/sub-' subj '_task-objects_run-*_events.tsv']);
    nr_runs = length(dataName);

    for data_nr = 1:nr_runs
        %%%% load preprocessed data
        load([dataRootPath '/sub-' subj '/ses-01/ieeg/' dataName(data_nr).name]);

        %%%% load stimulus information
        stim = readtable([dataRootPath '/sub-' subj '/ses-01/ieeg/' stimName(data_nr).name],...
            'FileType','text','Delimiter','\t');

        %%%% notch filter data at 60, 120 and 180 Hz
        data = ecog_notch(data,srate);

        %%%% make trial epochs
        onset_trial = round(stim.onset*srate);%from seconds to samples
        epoch_l = 1.5; % epoch length: -0.2:1.3 sec
        data_epoch = zeros(size(data,2),length(onset_trial),epoch_l*srate);
        for elec = 1:size(data,2)
            for l = 1:length(onset_trial)
                data_epoch(elec,l,:) = ...
                    data(onset_trial(l)-.2*srate+1:onset_trial(l)+(epoch_l-.2)*srate,elec)';
            end
        end
        % define t - time vector for each epoch
        t = [1:epoch_l*srate]/srate - 0.2;  
        
        %%%% make isi epochs
        onset_trial_off = round(stim.onset*srate)+1*srate;%from seconds to samples
        data_epoch_off = zeros(size(data,2),length(onset_trial_off),epoch_l*srate);
        for elec = 1:size(data,2)
            for l = 1:length(onset_trial_off)
                data_epoch_off(elec,l,:) = ...
                    data(onset_trial_off(l)-.2*srate+1:onset_trial_off(l)+(epoch_l-.2)*srate,elec)';
            end
        end

        %%%% baseline correct
%         data_epoch = ecog_baselinesubtract(data_epoch,t>-.1 & t<0);
%         data_epoch_off = ecog_baselinesubtract(data_epoch_off,t>-.1 & t<0);

        %%%% calculate spectra
        fft_w = window(@hann,100); % window width
        fft_ov = 50; % overlap
        % do not regress ERP here, because regressing out average response with a few trials can hurt 
        reg_erp = 0; % 1 to regress erp out, 0 not to, if yes, make sure to baseline correct first    

        % stimuli epochs
        fft_t = t>0 & t<=1; % time segment for spectrum % stimulus presentation 1 sec
        [f,data_epoch_spectra] = ...
            ecog_spectra(data_epoch,stim.trial_type,fft_w,fft_t,fft_ov,srate,reg_erp);
        stim_epoch = [stim.trial_type stim.stim_file_index];

        % ISI
        fft_t = t>0 & t<=.6; % time segment for spectrum % isi varied from 600-1400 ms
        [f,data_epoch_spectra_off] = ...
            ecog_spectra(data_epoch_off,ones(size(data_epoch_off,2),1),fft_w,fft_t,fft_ov,srate,reg_erp);

        disp(['done run ' int2str(data_nr)])

        % SETTINGS: fft_w = window(@hann,100), fft_ov = 50, reg_erp = 0, fft_t = t>=0 & t<1
        % no baseline subtract
        saveName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_run-' int2str(data_nr) '_spectra100.mat'];

        % SETTINGS: fft_w = window(@hann,300), fft_ov = 100, reg_erp = 0, fft_t = t>=0 & t<1
        % with baseline subtract
%         saveName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_run-' int2str(data_nr) '_spectra.mat'];

        % SETTINGS: fft_w = window(@hann,200), fft_ov = 100, reg_erp = 0, fft_t = t>=0 & t<1
        % with baseline subtract
%         saveName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_run-' int2str(data_nr) '_spectra200.mat'];

        % SETTINGS: fft_w = window(@hann,500), fft_ov = 250, reg_erp = 0, fft_t = t>=0 & t<1
        % with baseline subtract
%         saveName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_run-' int2str(data_nr) '_spectra500.mat'];
               
        % SETTINGS: fft_w = window(@hann,200), fft_ov = 100, reg_erp = 0, fft_t = t>=0 & t<.5
        % no baseline subtract
%         saveName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_run-' int2str(data_nr) '_spectraShort200.mat'];

        % SETTINGS: fft_w = window(@hann,500), fft_ov = 250, reg_erp = 0, fft_t = t>=.2 & t<1
        % no baseline subtract
%         saveName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_run-' int2str(data_nr) '_spectraLate500.mat'];

        % SETTINGS: fft_w = window(@hann,1000), fft_ov = 500, reg_erp = 0, fft_t = t>=0 & t<1
        % no baseline subtract, ISI window(@hann,500), fft_ov = 250
%         saveName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_run-' int2str(data_nr) '_spectra1000.mat'];

        save(saveName,'f','data_epoch_spectra','stim_epoch','data_epoch_spectra_off','exclude_channels','include_channels')   

    end % end run loop

end % end subject loop

%%
%% combine runs 
%%

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [9];

% analysisType = 'spectra';
% analysisType = 'spectraRERP';
% analysisType = 'spectra500';
% analysisType = 'spectraRERP500';
% analysisType = 'spectra200';
% analysisType = 'spectraRERP200';
% analysisType = 'spectraShort200';
% analysisType = 'spectraLate500';
% analysisType = 'spectra1000';
analysisType = 'spectra100';

for s = 1:length(subjects)
    % subject name
    if subjects(s)<10
        subj = ['0' int2str(subjects(s))];
    else
        subj = int2str(subjects(s));
    end

    dataName = dir([dataRootPath '/sub-' subj '/ses-01/ieeg/sub-' subj '_task-objects_run-*_ieeg_preproc.mat']);
    nr_runs = length(dataName); clear dataName
   
    % initialize empty variables
    spectra = [];  % temporary variable for epoched data
    spectra_off = [];  % temporary variable for epoched data
    stims = [];       % temporary variable for stimuli list
    runs = [];       % temporary variable for stimuli list

    for data_nr = 1:nr_runs

        dataName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_run-' int2str(data_nr) '_' analysisType '.mat'];

        load(dataName,'f','data_epoch_spectra','stim_epoch','data_epoch_spectra_off','exclude_channels','include_channels')   

        % concatonate with previous run
        spectra = cat(2,spectra,data_epoch_spectra);
        spectra_off = cat(2,spectra_off,data_epoch_spectra_off);
        stims = cat(1,stims,stim_epoch);
        runs = cat(1,runs,data_nr*ones(size(stim_epoch)));
    end

    saveName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '.mat'];

    save(saveName,'f','spectra','spectra_off','stims','runs','exclude_channels','include_channels')   

end

%% load all data together
%%
clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
subjects = [9];
s = 1;
% analysisType = 'spectra';
% analysisType = 'spectra500';
% analysisType = 'spectraRERP';
% analysisType = 'spectraRERP500';
% analysisType = 'spectra200';
analysisType = 'spectraLate500';

if subjects(s)<10
    subj = ['0' int2str(subjects(s))];
else
    subj = int2str(subjects(s));
end

dataName = [dataRootPath '/sub-' subj '/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '.mat'];
load(dataName)   

%% plot one channel to check:

% channels to plot (sub-09: 68):
electrodes=[68]; 

cm = jet(4);

for k=1:length(electrodes)
    figure('Color',[1 1 1],'Position',[0 0 300 300]),hold on
    elec=electrodes(k);
    
    data_fft=squeeze(spectra(elec,:,:));
    data_fft_off=squeeze(spectra_off(elec,:,:));
    
    title(['el ' int2str(elec)])
    
    stim_code=stims(:,1);
    stim_code(stim_code>=10 & stim_code<20)=1;
    stim_code(stim_code>=20 & stim_code<30)=2;
    stim_code(stim_code>=30 & stim_code<40)=3;
    stim_code(stim_code>=40 & stim_code<50)=4;

    % fft
    for s = 1:max(stim_code)
        plot(log10(f),log10(mean(data_fft(stim_code==s,:),1)),'-','Color',cm(s,:),'LineWidth',2)
    end
    plot(log10(f),log10(mean(data_fft_off,1)),'k','LineWidth',2)
    
    fill([log10(57) log10(64) log10(64) log10(57)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    fill([log10(117) log10(124) log10(124) log10(117)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    fill([log10(177) log10(184) log10(184) log10(177)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    
    xlim([log10(5) log10(200)])
    set(gca,'XTick',[log10(10) log10(25) log10(50) log10(100) log10(200)])
    set(gca,'XTickLabel',{'10','25','50','100','200'})
    
    set(gcf,'PaperPositionMode','auto')
%     print('-dpng','-r300',['./figures/spectra/' subj '_el'  int2str(electr)])
%     print('-depsc','-r300',['./figures/spectra/' subj '_el'  int2str(electr)])
%     close all
end


%% Plot one channel to check one stimulus from non-repeats:

% channels to plot (sub-09: 68):
electrodes = [68]; 
s = 50;

stims_plot = unique(stims(:,2));

for k = 1:length(electrodes)
    figure('Color',[1 1 1],'Position',[0 0 300 300]),hold on
    elec=electrodes(k);
    
    data_fft=squeeze(spectra(elec,stims(:,2)==stims_plot(s),:));
    data_fft_off=squeeze(spectra_off(elec,:,:));
    
    title(['el ' int2str(elec) ' stim ' int2str(stim_nr)])
    
    % fft
    plot(log10(f),log10(data_fft),'r-','LineWidth',2)
    plot(log10(f),log10(mean(data_fft_off,1)),'k','LineWidth',2)
    
    fill([log10(57) log10(64) log10(64) log10(57)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    fill([log10(117) log10(124) log10(124) log10(117)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    fill([log10(177) log10(184) log10(184) log10(177)],[2 2 -2 -2],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    
    xlim([log10(20) log10(200)])
    set(gca,'XTick',[log10(10) log10(25) log10(50) log10(100) log10(200)])
    set(gca,'XTickLabel',{'10','25','50','100','200'})
    
    set(gcf,'PaperPositionMode','auto')
%     print('-dpng','-r300',['./figures/spectra/' subj '_el'  int2str(electr)])
%     print('-depsc','-r300',['./figures/spectra/' subj '_el'  int2str(electr)])
%     close all
end

%% plot one channel to check repeats:

% channels to plot (sub-09: 66 67 68 (69)):
elec = [68]; 

stims_repeat = [];
stims_nonrepeat = [];
cntr_rep = 0;
cntr_nonrep = 0;
stimnr_list = unique(stims(:,2));
for k = 1:length(stimnr_list)
    stim_trials = length(find(stims(:,2)==stimnr_list(k)));
    if stim_trials == 1
        cntr_nonrep = cntr_nonrep+1;
        stims_nonrepeat(cntr_nonrep) = stimnr_list(k);
    elseif stim_trials > 1
        cntr_rep = cntr_rep+1;
        stims_repeat(cntr_rep) = stimnr_list(k);
    end
end

figure

% For repeated items:
for k = 1:length(stims_repeat)
    thisStims = find(stims(:,2)==stims_repeat(k));

    data_fft = squeeze(spectra(elec,thisStims,:));
    data_base = squeeze(mean(spectra_off(elec,:,:),2))'; 
    
    subplot(8,10,k),hold on
    
    plot(log10(f),log10(data_fft),'r')
    plot(log10(f),log10(data_base),'k')

    xlim([log10(20) log10(200)])
    ylim([-2.5 3.5])
    set(gca,'XTick',[log10(25) log10(50) log10(100) log10(200)])
    set(gca,'XTickLabel',[])
    set(gca,'YTick',[])
end
