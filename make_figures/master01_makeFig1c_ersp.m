% This script will generate the time-frequency plots for example electrodes
% in Figure 1c from: 
% Hermes D, Petridou N, Kay K, Winawer J. 2019 An image-computable model
% for the stimulus selectivity of gamma oscillations. bioRxiv doi:
% https://doi.org/10.1101/583567
%
% dhermes 2019 UMC Utrecht

% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath(genpath('/Users/dora/Documents/git/ecogBasicCode/'));
addpath(genpath('/Users/dora/Documents/m-files/knkutils'));
addpath(genpath('/Users/dora/Documents/m-files/Chronux/'))

%% load data and epoch

subjects = {'19','24','1001'};

s = 3;
% subject name
subj = subjects{s};
dataName = dir([dataDir '/sub-' subj '/ses-01/ieeg/sub-' subj '_ses-01_task-soc_run-*_ieeg_preproc.mat']);
stimName = dir([dataDir '/sub-' subj '/ses-01/ieeg/sub-' subj '_ses-01_task-soc_run-*_events.tsv']);
nr_runs = length(dataName);

data_epoch_all = [];
stim_all = [];

for data_nr = 1:nr_runs
    %%%% load preprocessed data
    load([dataDir '/sub-' subj '/ses-01/ieeg/' dataName(data_nr).name]);
    if size(data,1)>size(data,2)
        data = data';
    end
        
    %%%% load stimulus information
    stim = readtable([dataDir '/sub-' subj '/ses-01/ieeg/' stimName(data_nr).name],...
        'FileType','text','Delimiter','\t');

    %%%% notch filter data at 60, 120 and 180 Hz
    if isequal(subj,'19')
        data = ecog_notch(data',srate,60);
        data = data';
    elseif isequal(subj,'24')
        data = ecog_notch(data',srate,60);
        data = data';
    elseif isequal(subj,'1001')
        disp('no notch filter for subject 3, because not much line noise')
    end    
    
    %%%% make epochs
    onset_trial = round(stim.onset*srate);%from seconds to samples
    if srate==1000
        epoch_l = 1.2; % epoch length: -0.2:1 sec
        pre_stim = .2;
    elseif srate==2048
        epoch_l = 1.25; % epoch length: -0.2:1 sec
        pre_stim = .25;
    end
    data_epoch = zeros(size(data,1),length(onset_trial),epoch_l*srate);
    for elec = 1:size(data,1)
        for l = 1:length(onset_trial)
            data_epoch(elec,l,:) = ...
                data(elec,onset_trial(l)-pre_stim*srate+1:onset_trial(l)+(epoch_l-pre_stim)*srate);
        end
    end
    % define t - time vector for each epoch
    t = [1:epoch_l*srate]/srate - pre_stim;  
    
    data_epoch_all = cat(2,data_epoch_all,data_epoch);
    stim_all = cat(1,stim_all,stim.trial_type);

end

data_epoch = data_epoch_all;
% clear data_epoch_all data


%%
%% Plot ersp for an electrode and some stimuli
%%

elec = 50; 

conds_plot = {[45],[83],87};
% conds_plot = {39,40,41,42,43,44,45,46,87}; % orientation

cm1 = [repmat([0 0 0],100,1)];
cm1(1:40,1) = [0.7]';
cm1(1:40,2) = [0.7:-0.6/39:0.1]';
cm1(1:40,3) = [0.7:-0.7/39:0]';
cm1(40:100,1) = [0.7:(1-0.7)/60:1]';
cm1(40:100,2) = [0.1:.9/60:1]';

cm2 = [repmat([0 0 0],100,1)];
cm2(1:30,3) = [0.7]';
cm2(1:30,1) = [0.7:-0.7/29:0]';
cm2(1:30,2) = [0.7:-0.7/29:0]';
cm2(30:100,3) = [0.7:(1-0.7)/70:1]';
cm2(30:100,2) = [0:1/70:1]';
cm = [cm2(end:-1:1,:); cm1];

movingwin = [.200 .05];
params.pad = -1;
params.tapers = [3 5];
params.fpass = [0 200];
params.Fs = srate;
params.trialave = 0;

data_temp = squeeze(data_epoch(elec,:,:));
data_baseline = squeeze(data_epoch(elec,stim_all==87,:));

% calculate baseline
params.trialave = 1;
[S1b,t_tf_b,f_b] = mtspecgramc(data_baseline(:,t>.25 & t<.5)',movingwin,params);
S1b = mean(S1b,1);

params.trialave = 0;
figure('Color',[1 1 1],'Position',[0 0 350 250])
%%%%% all responses:
for kk = 1:length(conds_plot)
    data2use = squeeze(data_temp(ismember(stim_all,conds_plot{kk}),:))';

    % calculate spectgram
    [S1,t_tf,f] = mtspecgramc(data2use,movingwin,params);
    t_tf = t_tf+t(1);
    % normalize wrt baseline
    S1 = nanmean(S1,3);
    S1 = S1./repmat(S1b,[size(S1,1),1]);

    subplot(1,length(conds_plot),kk)
    imagesc(t_tf,f,log10(S1)',[-1.5 1.5])
    axis xy
    colormap(cm)
    set(gca,'XTick',[0 .5])
    title(['stim ' int2str(conds_plot{kk})])
end

colorbar

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',fullfile(dataDir,'derivatives','spectra','timefreq',...
%     ['ersp_sub-' subj '_el' int2str(elec) '_stim' int2str(conds_plot{1}) '_notch']))
% print('-depsc','-r300',fullfile(dataDir,'derivatives','spectra','timefreq',...
%     ['ersp_sub-' subj '_el' int2str(elec) '_stim' int2str(conds_plot{1}) '_notch']))

