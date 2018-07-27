clear all

addpath('~/Documents/git/ecogBasicCode/')
addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);
addpath(genpath('/Users/dorahermes-miller/Documents/m-files/knkutils'));

%% filter

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [19,23,24];

for s = 1:length(subjects)
% subject name
subj = subjects(s);
dataName = dir([dataRootPath '/sub-' int2str(subj) '/ieeg/sub-' int2str(subj) '_task-soc_run-*_ieeg_preproc.mat']);
stimName = dir([dataRootPath '/sub-' int2str(subj) '/ieeg/sub-' int2str(subj) '_task-soc_run-*_events.tsv']);
nr_runs = length(dataName);

for data_nr = 1:nr_runs
    %%%% load preprocessed data
    load([dataRootPath '/sub-' int2str(subj) '/ieeg/' dataName(data_nr).name]);

    %%%% notch filter data at 60, 120 and 180 Hz
    data = ecog_notch(data,1000);
    
    %filter
    hband_sig=zeros(size(data));
    for el=1:size(data,2)
        disp(['filter el ' int2str(el) ' of ' int2str(size(data,2))])
        hband=[101 200];
        hband_sig1=butterpass_eeglabdata(data(:,el),hband,srate);
        hband_sig(:,el)=hband_sig1;

        clear hband_sig1
    end
    clear data
    
    % SETTINGS01 filter from 101-200 Hz: 
    saveName = [dataRootPath '/sub-' int2str(subj) '/derivatives/ieeg/sub-' int2str(subj) '_task-soc_run-' int2str(data_nr) '_filter01.mat'];
    
    save(saveName,'hband_sig','exclude_channels','include_channels')   
    
end
end


