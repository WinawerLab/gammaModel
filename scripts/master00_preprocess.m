clear all

addpath(genpath('~/Documents/git/ecogBasicCode/'))
addpath('~/Documents/git/JSONio/')
addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);

% list of functions used:
% subj_visual_info.m
% ecog_notch.m
% ecog_baselinesubtract.m

%% task-soc downsample and CAR

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [19,23,24];

for s = 2%:length(subjects)
    % subject name
    subj = subjects(s);
    dataName = dir([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/sub-' int2str(subj) '_ses-01_task-soc_run-*_ieeg.mat']);
    nr_runs = length(dataName);
    dataJsonName = dir([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/sub-' int2str(subj) '_ses-01_task-soc_run-*_ieeg.json']);
    dataChannelsName = dir([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/sub-' int2str(subj) '_ses-01_task-soc_run-*_channels.tsv']);
    
    %%%% THIS LOOP IS FOR EPOCHING DATA ACROSS RUNS
    for data_nr = 1%:nr_runs 
        %%%% load raw data
        disp(['loading data set ' int2str(data_nr)])
        load([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/' dataName(data_nr).name]);
         
        %%%% get channel info 
        disp(['loading channels.tsv set ' int2str(data_nr)])
        channel_info = readtable([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/' dataChannelsName(data_nr).name],...
            'FileType','text','Delimiter','\t');
        % make a vector of channels to exclude from common average
        % referencing
        exclude_channels = zeros(length(channel_info.status),1);
        for k = 1:length(channel_info.status)
            if isequal(channel_info.status{k},'bad')  
                exclude_channels(k) = 1;
            end
            if ~isequal(channel_info.type{k},'iEEG')  
                exclude_channels(k) = 1;
            end
        end
        include_channels = find(exclude_channels==0);
        exclude_channels = find(exclude_channels>0);
        
        %%%% get iEEG info to get sampling frequency
        disp(['loading ieeg.json set ' int2str(data_nr)])
        ieeg_json = jsonread([dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/' dataJsonName(data_nr).name]);
        srate = round(ieeg_json.SamplingFrequency);
        
        % plot a quick figure of data quality here
        figure
        % plot downsampled timeseries included and excluded channels
        % resample to 1 value/second - this is only for a quick dataview
        temp_data = resample(data,1000,1031);
        temp_data = resample(temp_data,1,1480);
        subplot(3,2,1),hold on
        title('quick plot of downsampled data, do not zoom')
        plot(temp_data(:,include_channels))
        subplot(3,2,3),hold on
        plot(temp_data(:,exclude_channels))
        clear temp_data
        % plot power spectra of included and excluded channels
        [pxx,f] = pwelch(data,2*window(@hann,srate),0,2*srate,srate);
        subplot(1,4,3),hold on
        plot(f,log10(pxx(:,include_channels)),'b')
        plot(f,log10(pxx(:,exclude_channels)),'r')
        xlim([0 300])
        
        % rereference
        disp(['CAR'])
        [data] = ecog_CarRegress(data,include_channels);
        
        % resample / cut down to 1000 Hz
        s = srate*(1000/1031)*(1000/1480); % check that it's appropriate    
        disp(['downsampling to ' int2str(s) ' Hz'])
        data = resample(data,1000,1031);
        data = resample(data,1000,1480);
        srate = floor(s); clear s
        disp('data downsampled')
        
        temp_data = resample(data,1,1000);
        subplot(3,2,5),hold on
        plot(temp_data(:,include_channels))
        clear temp_data
        
        % plot power spectra of included and excluded channels
        [pxx,f] = pwelch(data,2*window(@hann,srate),0,2*srate,srate);
        subplot(1,4,4),hold on
        plot(f,log10(pxx(:,include_channels)),'b')
        xlim([0 300])
        
        [~,a,~] = fileparts(dataName(data_nr).name);
        dataNamePrep = [dataRootPath '/sub-' int2str(subj) '/ses-01/ieeg/' a '_preproc.mat'];
        %%%% SAVE DATA
        save(dataNamePrep,'-v7.3','data','srate','include_channels','exclude_channels')
        clear data
    end

end


%% task-objects downsample and CAR

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

    dataName = dir([dataRootPath '/sub-' subj '/ieeg/sub-' subj '_task-objects_run-*_ieeg.mat']);
    nr_runs = length(dataName);
    dataJsonName = dir([dataRootPath '/sub-' subj '/ieeg/sub-' subj '_task-objects_run-*_ieeg.json']);
    dataChannelsName = dir([dataRootPath '/sub-' subj '/ieeg/sub-' subj '_task-objects_run-*_channels.tsv']);
    
    %%%% THIS LOOP IS FOR EPOCHING DATA ACROSS RUNS
    for data_nr = 1:nr_runs 
        %%%% load raw data
        disp(['loading data set ' int2str(data_nr)])
        load([dataRootPath '/sub-' subj '/ieeg/' dataName(data_nr).name]);
         
        %%%% get channel info 
        disp(['loading channels.tsv set ' int2str(data_nr)])
        channel_info = readtable([dataRootPath '/sub-' subj '/ieeg/' dataChannelsName(data_nr).name],...
            'FileType','text','Delimiter','\t');
        % make a vector of channels to exclude from common average
        % referencing
        exclude_channels = zeros(length(channel_info.status),1);
        for k = 1:length(channel_info.status)
            if isequal(channel_info.status{k},'bad')  
                exclude_channels(k) = 1;
            end
            if ~isequal(channel_info.type{k},'iEEG')  
                exclude_channels(k) = 1;
            end
        end
        include_channels = find(exclude_channels==0);
        exclude_channels = find(exclude_channels>0);
        
        %%%% get iEEG info to get sampling frequency
        disp(['loading ieeg.json set ' int2str(data_nr)])
        ieeg_json = jsonread([dataRootPath '/sub-' subj '/ieeg/' dataJsonName(data_nr).name]);
        srate = round(ieeg_json.SamplingFrequency);
        
        % plot a quick figure of data quality here
        figure
        % plot downsampled timeseries included and excluded channels
        % resample to 1 value/second - this is only for a quick dataview
        temp_data = resample(data,1000,1031);
        temp_data = resample(temp_data,1,2960);
        subplot(3,2,1),hold on
        title('quick plot of downsampled data, do not zoom')
        plot(temp_data(:,include_channels))
        subplot(3,2,3),hold on
        plot(temp_data(:,exclude_channels))
        clear temp_data
        % plot power spectra of included and excluded channels
        [pxx,f] = pwelch(data,2*window(@hann,srate),0,2*srate,srate);
        subplot(1,4,3),hold on
        plot(f,log10(pxx(:,include_channels)),'b')
        plot(f,log10(pxx(:,exclude_channels)),'r')
        xlim([0 300])
        xlabel('Frequency (Hz)')
        title('raw data')
        
        % rereference
        disp(['CAR'])
        [data] = ecog_CarRegress(data,include_channels);
        
        % resample / cut down to 1000 Hz
        s=srate*(1000/1031)*(1000/2960); 
        disp(['downsampling to ' int2str(s) ' Hz'])
        data=resample(data,1000,1031);
        data=resample(data,1000,2960);
        srate=floor(s); clear s
        disp('data downsampled')
        
        temp_data = resample(data,1,1000);
        subplot(3,2,5),hold on
        plot(temp_data(:,include_channels))
        clear temp_data
        
        % plot power spectra of included and excluded channels
        [pxx,f] = pwelch(data,2*window(@hann,srate),0,2*srate,srate);
        subplot(1,4,4),hold on
        plot(f,log10(pxx(:,include_channels)),'b')
        xlim([0 300])
        xlabel('Frequency (Hz)')
        title('good channels after CAR and downsample')
        
        [~,a,~] = fileparts(dataName(data_nr).name);
        dataNamePrep = [dataRootPath '/sub-' subj '/ieeg/' a '_preproc.mat'];
        %%%% SAVE DATA
        save(dataNamePrep,'-v7.3','data','srate','include_channels','exclude_channels')
        clear data
    end

end
