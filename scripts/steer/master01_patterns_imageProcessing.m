
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

% load images
load(fullfile(dataRootPath,'stimuli','task-soc_stimuli.mat'),'stimuli')

% downsample to increase speed
temp = zeros(240,240,size(stimuli,3),'single');
for p=1:size(stimuli,3)
  temp(:,:,p) = imresize(single(stimuli(:,:,p)),[240 240],'cubic');
end
stimulus = temp;
clear temp;

% im = double(stimulus(:,:,10));
im = double(stimulus(:,:,11));

% Divide by mean and subtract 1 so that it is a "contrast"
% image.
im = im/mean2(im) - 1;

%% Notes on spatial frequency:

% 800x800 images result in numLevels = 7
% Level 7: 1 cycle per image: 1 cycle ~ 20 degrees
% Level 6: 2 cycles per image: 1 cycle ~ 10 degrees
% Level 5: 4 cycles per image: 1 cycle ~ 5 degrees
% Level 4: 8 cycles per image: 1 cycle ~ 2.5 degrees
% Level 3: 16 cycles per image: 1 cycle ~ 1.25 degrees
% Level 2: 32 cycles per image: 1 cycle ~ 0.625 degrees
% Level 1: 64 cycles per image: 1 cycle ~ 0.31 degrees

% Images downsampled to 240 x 240 result in numLevels = 5
% Level 5: 1 cycle per image: 1 cycle ~ 20 degrees
% Level 4: 2 cycles per image: 1 cycle ~ 10 degrees
% Level 3: 4 cycles per image: 1 cycle ~ 5 degrees
% Level 2: 8 cycles per image: 1 cycle ~ 2.5 degrees
% Level 1: 16 cycles per image: 1 cycle ~ 1.25 degrees

%% Filter one image to test
tic
% Settings for steerable pyramid filters
numOrientations = 8;
bandwidth       = 1; 
dims            = size(im); 
numLevels       = maxLevel(dims,bandwidth); % calculates nr of spatial frequencies

% Design filters
[freqRespsImag,freqRespsReal,pind]= ...
    makeQuadFRs(dims,numLevels,numOrientations,bandwidth); % pind will be size: (numLevels x numOrientations) + 2
plot_orient = 1;
plot_level = 1;
a = accessSteerBand(freqRespsImag,pind,numOrientations,plot_level,plot_orient);
% figure,imagesc(abs(a))
for plot_orient = 1:numOrientations
    a = accessSteerBand(freqRespsImag,pind,numOrientations,plot_level,plot_orient);
    figure,displayImage(a);
%     set(gcf,'PaperPositionMode','auto')
%     print('-depsc','-r300',['./figures/steerFilters/ExampleFilters_Level' int2str(plot_level) '_Or' int2str(plot_orient)])
%     print('-dpng','-r300',['./figures/steerFilters/ExampleFilters_Level' int2str(plot_level)  '_Or' int2str(plot_orient)])
end

% visualize filters
% viewBands(freqRespsImag,pind,1/4,'auto1');
% viewBands(freqRespsReal,pind,1/4,'auto1');

% Run filters through image
[pyr,pind]=buildQuadBands(im,freqRespsImag,freqRespsReal);
clear freqRespsImag freqRespsReal % housekeeping
% visualine pyr (pyr is complex, visualize real):
plot_orient = 1;
plot_level = 2;
displayImage(accessSteerBand(pyr,pind,numOrientations,plot_level,plot_orient));
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

% Save weighted sum with weights centered at one of 3 levels:
test_out1 = zeros(dims(1),dims(2),numOrientations); % 2: high sf
test_out2 = zeros(dims(1),dims(2),numOrientations); % 4: mid sf
test_out3 = zeros(dims(1),dims(2),numOrientations); % 6: low sf

% Get weights:
w = window(@gausswin,5);
w1 = [w(2:5); zeros(3,1)];
w2 = [0; w; 0];
w3 = [zeros(3,1); w(1:4)];
for plot_orient = 1:numOrientations
    temp_band = zeros(numLevels,dims(1),dims(2));
    for plot_level = 1:numLevels
        temp_band(plot_level,:,:) = accessSteerBand(nEnergies,pind,numOrientations,plot_level,plot_orient);
    end
    test_out1(:,:,plot_orient) = squeeze(sum(temp_band*w1,1));
    test_out2(:,:,plot_orient) = squeeze(sum(temp_band*w2,1));
    test_out3(:,:,plot_orient) = squeeze(sum(temp_band*w3,1));
    clear temp_band
end

toc

%% Loop through images and save
clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

% load images
load(fullfile(dataRootPath,'stimuli','task-soc_stimuli.mat'),'stimuli')

% downsample to increase speed
temp = zeros(240,240,size(stimuli,3),'single');
for p=1:size(stimuli,3)
  temp(:,:,p) = imresize(single(stimuli(:,:,p)),[240 240],'cubic');
end
stimuli = temp;
clear temp;

% settings for steerable pyramid filters
numOrientations = 8;
bandwidth       = 1; 
dims            = [size(stimuli,1), size(stimuli,2)]; 
numLevels       = maxLevel(dims,bandwidth); % calculates nr of spatial frequencies

% predefine output
imEnergy = single(zeros(size(stimuli,3),size(stimuli,1)*size(stimuli,2)*numOrientations,numLevels));

for kk = 1:size(stimuli,3)
    im = double(stimuli(:,:,kk));

    % Divide by mean and subtract 1 so that it is a "contrast"
    % image.
    im = im/mean2(im) - 1;   
    
    % Design filters
    [freqRespsImag,freqRespsReal,pind] = ...
        makeQuadFRs(dims,numLevels,numOrientations,bandwidth); % pind will be size: (numLevels x numOrientations) + 2

    % Run filters through image
    [pyr,pind]=buildQuadBands(im,freqRespsImag,freqRespsReal);
    clear freqRespsImag freqRespsReal % housekeeping

    % Calculate energy from filtered images (pyr):
    nEnergies = normEnergies(pyr,pind,numOrientations,0.1);

    % Put spatial frequencies/orientations in output imEnergy 
    for plot_level = 1:numLevels
        temp_band = zeros(numOrientations,dims(1),dims(2));
        for plot_orient = 1:numOrientations
            temp_band(plot_orient,:,:) = accessSteerBand(nEnergies,pind,numOrientations,plot_level,plot_orient);
        end
        imEnergy(kk,:,plot_level) = temp_band(:);
    end 
    
end

% save image energies
save(['./data/task-soc_stimuli_steerFilt01.mat'],'imEnergy','-v7.3' )

%%
% load(['./data/task-soc_stimuli_steerFilt01.mat'],'imEnergy')
% numOrientations = 8;
res = sqrt(size(imEnergy,2)./numOrientations);

% Now to view an image, one can use:
inNr = 10;
plot_level = 1;
thisIm = reshape(imEnergy(inNr,:,plot_level),numOrientations,res,res); 
figure('Position',[0 0 1200 300]), 
for kk = 1:numOrientations
    subplot(1,numOrientations+1,kk)
    imagesc(squeeze(thisIm(kk,:,:)),[0 max(thisIm(:))/2])
    axis square
    axis off
end
colormap gray

inNr = 10;
plot_level = 1;
thisIm = reshape(imEnergy(inNr,:,plot_level),numOrientations,res,res); 
subplot(1,numOrientations+1,numOrientations+1),imagesc(squeeze(mean(thisIm,1)),[0 max(thisIm(:))/2])
axis square
axis off
title('mean across orientations')

set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/steerFilters/ExampleStim' int2str(inNr) '_Level' int2str(plot_level)])
% print('-dpng','-r300',['./figures/steerFilters/ExampleStim' int2str(inNr) '_Level' int2str(plot_level)])

%% plot one image to check size

load(['./data/task-soc_stimuli_steerFilt01.mat'],'imEnergy')
numOrientations = 8;
res = sqrt(size(imEnergy,2)/numOrientations);

% Now to view an image, one can use:
inNr = 11;
plot_level = 1;

figure('Position',[0 0 1200 200]) 

thisIm = reshape(imEnergy(inNr,:,plot_level),numOrientations,res,res); 

for kk = 1:numOrientations
    subplot(1,numOrientations+1,kk)
    imagesc(squeeze(thisIm(kk,:,:)),[0 max(thisIm(:))/2])
    axis square
    axis off
end
colormap gray

thisIm = reshape(imEnergy(inNr,:,plot_level),numOrientations,res,res); 
subplot(1,numOrientations+1,numOrientations+1),imagesc(squeeze(mean(thisIm,1)),[0 .1])
colorbar
axis square
axis off
title('mean across orientations')

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/steerFilters/ExampleStim_ds_' int2str(inNr) '_Level' int2str(plot_level)])
% print('-dpng','-r300',['./figures/steerFilters/ExampleStim_ds_' int2str(inNr) '_Level' int2str(plot_level)])

