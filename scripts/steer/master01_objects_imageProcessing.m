
%
% Filter patterns with steerable pyramid filters. Original code was
% developed by David Heeger, here adapted to filter patterns from ECoG
% experiments.
%
% Dora Hermes, 2017

addpath('~/Documents/git/ecogBasicCode/')
addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);
addpath(genpath('~/Documents/m-files/knkutils'));
addpath('/Users/dora/Documents/m-files/steer/')

%%
%% %%%%%% LOAD images %%%%%%%%
%%

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

% load images - these are resampled as in
% /m-files/model/master01_objects_imageProcessing.m

subj = '09';
load(['../model/data/sub-' subj '_stimuli'],'ims_repeat','ims_nonrepeat')

% Divide by mean and subtract 1 so that it is a "contrast"
% image.
im = ims_repeat(:,:,1);
im(im<0)=0; % already the case, just to make sure
im(im>254)=254; % already the case, just to make sure
im = im/mean2(im) - 1; % scale  


%% Example for first image to test code
tic
% Settings for steerable pyramid filters
numOrientations = 8;
bandwidth       = 1; 
dims            = size(im); 
numLevels       = maxLevel(dims,bandwidth); % calculates nr of spatial frequencies

% Design filters
[freqRespsImag,freqRespsReal,pind]= ...
    makeQuadFRs(dims,numLevels,numOrientations,bandwidth); % pind will be size: (numLevels x numOrientations) + 2
% visualize filters
% viewBands(freqRespsImag,pind,1/4,'auto1');
% viewBands(freqRespsReal,pind,1/4,'auto1');

% Run filters through image
[pyr,pind]=buildQuadBands(im,freqRespsImag,freqRespsReal);
clear freqRespsImag freqRespsReal % housekeeping
% visualine pyr:
% plot_orient = 1;
% plot_level = 7;
% displayImage(accessSteerBand(pyr,pind,numOrientations,plot_level,plot_orient));
% displayImage(pyrLow(pyr,pind));
% displayImage(pyrHi(pyr,pind));
% viewBands(pyr,pind,1/4,'auto1');

% Calculate energy from filtered images (pyr):
nEnergies = normEnergies(pyr,pind,numOrientations,0.1);
max2(nEnergies)
% visualize energies:
% viewBands(nEnergies,pind,1/4,[0 0.8]); % max possible response is 0.8
% viewBands(nEnergies,pind,1/4,[0 max2(nEnergies)]);
% plot_orient = 3;
% plot_level = 3;
% band = accessSteerBand(nEnergies,pind,numOrientations,plot_level,plot_orient);
% displayImage(band,[0 0.8]);
% imStats(band);
%%
% Just sum across spatial frequencies? 
test_out = zeros(dims(1),dims(2),numOrientations);
for plot_orient = 1:numOrientations
    temp_band = zeros(numLevels,dims(1),dims(2));
    for plot_level = 1:numLevels
        temp_band(plot_level,:,:) = accessSteerBand(nEnergies,pind,numOrientations,plot_level,plot_orient);
    end
    test_out(:,:,plot_orient) = squeeze(sum(temp_band,1));
    clear temp_band
end

toc

%% Loop through repeated images and save
clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

% load images
subj = '09';
load(['../model/data/sub-' subj '_stimuli'],'ims_repeat','ims_nonrepeat')

stimuli = ims_repeat;

% settings
numOrientations = 4;
bandwidth       = 1; 

% predefine output
imEnergy = single(zeros(size(stimuli,3),size(stimuli,1)*size(stimuli,2)*numOrientations));
imEnergyHigh = single(zeros(size(stimuli,3),size(stimuli,1)*size(stimuli,2)*numOrientations));
imEnergyMid = single(zeros(size(stimuli,3),size(stimuli,1)*size(stimuli,2)*numOrientations));
imEnergyLow = single(zeros(size(stimuli,3),size(stimuli,1)*size(stimuli,2)*numOrientations));

for kk = 1:size(stimuli,3)
    im = double(stimuli(:,:,kk));

    % Divide by mean and subtract 1 so that it is a "contrast"
    % image.
    im = im/mean2(im) - 1;
%     im = im/mean(stimuli(:)) - 1;
    
    % other settings for steerable pyramid filters
    dims            = size(im); 
    numLevels       = maxLevel(dims,bandwidth); % calculates nr of spatial frequencies
    
    % Design filters
    [freqRespsImag,freqRespsReal,pind] = ...
        makeQuadFRs(dims,numLevels,numOrientations,bandwidth); % pind will be size: (numLevels x numOrientations) + 2

    % Run filters through image
    [pyr,pind]=buildQuadBands(im,freqRespsImag,freqRespsReal);
    clear freqRespsImag freqRespsReal % housekeeping

    % Calculate energy from filtered images (pyr):
    nEnergies = normEnergies(pyr,pind,numOrientations,0.1);

    % Sum across spatial frequencies (not sure if this is the right thing to do) 
    test_out = zeros(dims(1),dims(2),numOrientations);
    for plot_orient = 1:numOrientations
        temp_band = zeros(numLevels,dims(1),dims(2));
        for plot_level = 1:numLevels
            temp_band(plot_level,:,:) = accessSteerBand(nEnergies,pind,numOrientations,plot_level,plot_orient);
        end
        test_out(:,:,plot_orient) = squeeze(sum(temp_band,1));
        clear temp_band
    end
    imEnergy(kk,:) = single(test_out(:));
    
    % Save weighted sum with weights centered at one of 3 levels:
    test_out1 = zeros(dims(1),dims(2),numOrientations); % 2: high sf
    test_out2 = zeros(dims(1),dims(2),numOrientations); % 4: mid sf
    test_out3 = zeros(dims(1),dims(2),numOrientations); % 6: low sf

    % Get weights:
    w = window(@gausswin,5);
    w1 = [w(3:5); zeros(3,1)];
    w2 = [w; 0];
    w3 = [zeros(2,1); w(1:4)];
    for plot_orient = 1:numOrientations
        temp_band = zeros(numLevels,dims(1),dims(2));
        for plot_level = 1:numLevels
            temp_band(plot_level,:,:) = accessSteerBand(nEnergies,pind,numOrientations,plot_level,plot_orient);
        end
        test_out1(:,:,plot_orient) = squeeze(sum(temp_band.*w1,1));
        test_out2(:,:,plot_orient) = squeeze(sum(temp_band.*w2,1));
        test_out3(:,:,plot_orient) = squeeze(sum(temp_band.*w3,1));
        clear temp_band
    end
    imEnergyHigh(kk,:) = single(test_out1(:));
    imEnergyMid(kk,:) = single(test_out2(:));
    imEnergyLow(kk,:) = single(test_out3(:));
end

% save image energies
% 01: Using: im = im/mean2(im) - 1;
save(['./data/task-objects_stimuli-repeat_steerFilt01.mat'],'imEnergy')
save(['./data/task-objects_stimuli-repeat_steerFilt01_high.mat'],'imEnergyHigh')
save(['./data/task-objects_stimuli-repeat_steerFilt01_mid.mat'],'imEnergyMid')
save(['./data/task-objects_stimuli-repeat_steerFilt01_low.mat'],'imEnergyLow')

% 02: Using: im = im/mean(stimuli(:)) - 1;
% save(['./data/task-objects_stimuli-repeat_steerFilt02.mat'],'imEnergy')

% now to view an image, one can use:
im_nr = 1;
a = reshape(imEnergy(im_nr,:),dims(1),dims(2),4);
figure,imagesc(a(:,:,4))

%% Loop through nonrepeated images and save
clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

% load images
subj = '09';
load(['../model/data/sub-' subj '_stimuli'],'ims_repeat','ims_nonrepeat')

stimuli = ims_nonrepeat;

% settings
numOrientations = 4;
bandwidth       = 1; 

% predefine output
imEnergy = single(zeros(size(stimuli,3),size(stimuli,1)*size(stimuli,2)*numOrientations));
imEnergyHigh = single(zeros(size(stimuli,3),size(stimuli,1)*size(stimuli,2)*numOrientations));
imEnergyMid = single(zeros(size(stimuli,3),size(stimuli,1)*size(stimuli,2)*numOrientations));
imEnergyLow = single(zeros(size(stimuli,3),size(stimuli,1)*size(stimuli,2)*numOrientations));

for kk = 1:size(stimuli,3)
    im = double(stimuli(:,:,kk));

    % Divide by mean and subtract 1 so that it is a "contrast"
    % image.
    im = im/mean2(im) - 1;
    
    % other settings for steerable pyramid filters
    dims            = size(im); 
    numLevels       = maxLevel(dims,bandwidth); % calculates nr of spatial frequencies
    
    % Design filters
    [freqRespsImag,freqRespsReal,pind] = ...
        makeQuadFRs(dims,numLevels,numOrientations,bandwidth); % pind will be size: (numLevels x numOrientations) + 2

    % Run filters through image
    [pyr,pind]=buildQuadBands(im,freqRespsImag,freqRespsReal);
    clear freqRespsImag freqRespsReal % housekeeping

    % Calculate energy from filtered images (pyr):
    nEnergies = normEnergies(pyr,pind,numOrientations,0.1);

    % Sum across spatial frequencies (not sure if this is the right thing to do) 
    test_out = zeros(dims(1),dims(2),numOrientations);
    for plot_orient = 1:numOrientations
        temp_band = zeros(numLevels,dims(1),dims(2));
        for plot_level = 1:numLevels
            temp_band(plot_level,:,:) = accessSteerBand(nEnergies,pind,numOrientations,plot_level,plot_orient);
        end
        test_out(:,:,plot_orient) = squeeze(sum(temp_band,1));
        clear temp_band
    end
    imEnergy(kk,:) = single(test_out(:));
    
     % Save weighted sum with weights centered at one of 3 levels:
    test_out1 = zeros(dims(1),dims(2),numOrientations); % 2: high sf
    test_out2 = zeros(dims(1),dims(2),numOrientations); % 4: mid sf
    test_out3 = zeros(dims(1),dims(2),numOrientations); % 6: low sf

    % Get weights:
    w = window(@gausswin,5);
    w1 = [w(3:5); zeros(3,1)];
    w2 = [w; 0];
    w3 = [zeros(2,1); w(1:4)];
    for plot_orient = 1:numOrientations
        temp_band = zeros(numLevels,dims(1),dims(2));
        for plot_level = 1:numLevels
            temp_band(plot_level,:,:) = accessSteerBand(nEnergies,pind,numOrientations,plot_level,plot_orient);
        end
        test_out1(:,:,plot_orient) = squeeze(sum(temp_band.*w1,1));
        test_out2(:,:,plot_orient) = squeeze(sum(temp_band.*w2,1));
        test_out3(:,:,plot_orient) = squeeze(sum(temp_band.*w3,1));
        clear temp_band
    end
    imEnergyHigh(kk,:) = single(test_out1(:));
    imEnergyMid(kk,:) = single(test_out2(:));
    imEnergyLow(kk,:) = single(test_out3(:));
end

% save image energies
save(['./data/task-objects_stimuli-nonrepeat_steerFilt01.mat'],'imEnergy')
save(['./data/task-objects_stimuli-nonrepeat_steerFilt01_high.mat'],'imEnergyHigh')
save(['./data/task-objects_stimuli-nonrepeat_steerFilt01_mid.mat'],'imEnergyMid')
save(['./data/task-objects_stimuli-nonrepeat_steerFilt01_low.mat'],'imEnergyLow')

% now to view an image, one can use:
im_nr = 1;
a = reshape(imEnergy(im_nr,:),dims(1),dims(2),4);
figure,imagesc(a(:,:,4))