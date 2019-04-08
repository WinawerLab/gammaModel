clear all

% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath(genpath('/Users/dora/Documents/git/ecogBasicCode/'));
addpath(genpath('/Users/dora/Documents/m-files/knkutils'));
addpath('~/Documents/git/JSONio/')
% also use Fieldtrip to read in the data

% list of functions used:
% subj_visual_info.m
% ecog_notch.m
% ecog_baselinesubtract.m

%% task-soc downsample and CAR

subjects = {'19','24','1001'};

for s = 1%:length(subjects)
    % subject name
    subj = subjects{s};
    dataName = dir(fullfile(dataDir,['sub-' subj],['ses-01'],'ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_ieeg.eeg']));
    nr_runs = length(dataName);
    dataJsonName = dir(fullfile(dataDir,['sub-' subj],['ses-01'],'ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_ieeg.json']));
    dataChannelsName = dir(fullfile(dataDir,['sub-' subj],['ses-01'],'ieeg',...
        ['sub-' subj '_ses-01_task-soc_run-*_channels.tsv']));
    
    %%%% THIS LOOP IS FOR PREPROCESSING DATA ACROSS RUNS
    for data_nr = 1:nr_runs 
        %%%% load raw data
        disp(['loading data set ' int2str(data_nr)])
        data = ft_read_data(fullfile(dataDir,['sub-' subj],'ses-01','ieeg',...
            dataName(data_nr).name));
        data_hdr = ft_read_header(fullfile(dataDir,['sub-' subj],'ses-01','ieeg',...
            dataName(data_nr).name));
         srate = round(data_hdr.Fs); % get sampling frequency

        %%%% get channel info 
        disp(['loading channels.tsv set ' int2str(data_nr)])
        channel_info = readtable(fullfile(dataDir,['sub-' subj],'ses-01','ieeg/',...
            dataChannelsName(data_nr).name),...
            'FileType','text','Delimiter','\t');
        % make a vector of channels to exclude from common average
        % referencing
        exclude_channels = zeros(length(channel_info.status),1);
        for k = 1:length(channel_info.status)
            if isequal(channel_info.status{k},'bad')  
                exclude_channels(k) = 1;
            end
            if ~isequal(channel_info.type{k},'ECOG') && ~isequal(channel_info.type{k},'SEEG')
                exclude_channels(k) = 1;
            end
        end
        include_channels = find(exclude_channels==0);
        exclude_channels = find(exclude_channels>0);
        
        % rereference
        disp(['CAR'])
        [data] = ecog_CarRegress(data,include_channels);
       
        % resample / cut awkward sampling rates down to 1000 Hz
        if ismember('24',{'19','24'})
            s = srate*(1000/1031)*(1000/1480); 
            disp(['downsampling to ' int2str(s) ' Hz'])
            data = resample(data,1000,1031);
            data = resample(data,1000,1480);
            srate = floor(s); clear s
            disp('data downsampled')
        end
        
        [~,a,~] = fileparts(dataName(data_nr).name);
        dataNamePrep = fullfile(dataDir,'derivatives','preprocessing',['sub-' subj],'ses-01','ieeg',...
            [a '_preproc.mat']);
        %%%% SAVE DATA
        save(dataNamePrep,'-v7.3','data','srate','include_channels','exclude_channels')
        clear data
    end

end

