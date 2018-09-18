% Stand alone description of the NBF model
%
% Processing steps
% 1. Images are filtered with Gabor patches
% 2. Squaring and summing across phases to get complex cell energy.
% 3. Spatial summation across a population receptive field.
% 4. Standard deviation across orientation.
%
% The original code to get complex cell energy was developed by Kendrick
% Kay (PLOS Computational Biology, 2013).
%
% Dora Hermes, UMC Utrecht, 2018.

%%
% Adding toolbox for filtering the images with Gabor patches:
% from https://github.com/kendrickkay/knkutils
addpath(genpath('~/Documents/m-files/knkutils'));

% Load images (these are stored in a .mat file)
load('/Volumes/DoraBigDrive/data/visual_soc/soc_bids/stimuli/task-soc_stimuli.mat','stimuli')

%% Downsample and zero-pad images

%%%% DOWNSAMPLE TO INCREASE SPEED
% resize the stimuli to 240 x 240 to reduce computational time.
% use single-format to save memory.
temp = zeros(240,240,size(stimuli,3),'single');
for p=1:size(stimuli,3)
  temp(:,:,p) = imresize(single(stimuli(:,:,p)),[240 240],'cubic');
end
stimulus = temp;
clear temp;

%%%% RESCALE
% ensure that all values are between 0 and 254.
% rescale values to the range [0,1].
% subtract off the background luminance (0.5).
% after these steps, all pixel values will be in the
% range [-.5,.5] with the background corresponding to 0.
stimulus(stimulus < 0) = 0;
stimulus(stimulus > 254) = 254;
stimulus = stimulus/254 - 0.5;

%%%% ZERO PAD
% pad the stimulus with zeros to reduce edge effects.
% the new resolution is 270 x 250 (15-pixel padding on each side).
stimulus = placematrix(zeros(270,270,size(stimulus,3),'single'),stimulus);

% inspect one of the stimuli
figure;
imagesc(stimulus(:,:,15));
axis image tight;
caxis([-.5 .5]);
colormap(gray);
colorbar;
title('Stimulus');

%% 1. Filtering images with Gabor patches 
% takes ~10 mins for 86 images %

% Apply Gabor filters to the stimuli. Filters occur at different positions,
% orientations, and phases. There are several parameters that govern the
% design of the filters:
filt_prop.cycles = 60*(270/240);    %   the number of cycles per image is 60*(270/240)
filt_prop.bandwidth = -1;           %   the spatial frequency bandwidth of the filters is 1 octave
filt_prop.spacings = 1;             %   the separation of adjacent filters is 1 std dev of the Gaussian envelopes
                                    %     (this results in a 90 x 90 grid of positions)
filt_prop.orientations = 8;         %   filters occur at 8 orientations
filt_prop.phases = 2;               %   filters occur at 2 phases (between 0 and pi)
filt_prop.thres = 0.01;             %   the Gaussian envelopes are thresholded at .01
filt_prop.scaling = 2;              %   filters are scaled to have an equivalent Michelson contrast of 1
filt_prop.mode = 0;                 %   the dot-product between each filter and each stimulus is computed

% After this step, stimulus is images x phases*orientations*positions.
stimuliFilt = applymultiscalegaborfilters(reshape(stimulus,270*270,[])', ...
  filt_prop.cycles,filt_prop.bandwidth,filt_prop.spacings,filt_prop.orientations,...
  filt_prop.phases,filt_prop.thres,filt_prop.scaling,filt_prop.mode);

%% 2. Squaring and summing across phases to get complex cell energy.

stimuliEnergy = sqrt(blob(stimuliFilt.^2,2,2)); % quadrature pairs

% Now to make a figure of the energy for each orientation for one image,
% one can use:
figure('Position',[0 0 1200 300]), 
res = sqrt(size(stimuliEnergy,2)/filt_prop.orientations);
inNr = 12;
% reshape to orientationXiXj
thisIm = reshape(stimuliEnergy(inNr,:),filt_prop.orientations,res,res); 
% plot each orientation in a subplot
for kk = 1:numOrientations
    subplot(1,numOrientations+1,kk)
    imagesc(squeeze(thisIm(kk,:,:)),[0 max(thisIm(:))])
    axis square
    axis off
end
colormap gray
% Add lastplot of the mean across orientations
subplot(1,numOrientations+1,numOrientations+1),imagesc(squeeze(mean(thisIm,1)),[0 max(thisIm(:))/2])
axis square
axis off
title('mean across orientations')

%% 3. Spatial summation across a certain area in the image.

% Get resolution of the stimulus:
res = sqrt(size(stimuliEnergy,2)/filt_prop.orientations);

% Let's sum across a center area, this sets the center and size of a
% Gaussian across which to sum:
pp = [res/2 res/2 res/10];

% Issue a dummy call to makegaussian2d.m to pre-compute xx and yy.
% these variables are re-used to achieve faster computation.
[~,xx,yy] = makegaussian2d(res,2,2,2,2);

% Define a function that computes a 2D Gaussian across which to sum. The
% Gaussian has center pp(1), pp(2) and size pp(3). The ouput is a vector of
% size res*res x 1
gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

% Sum across the Gaussian:
stimuliEnergyPRF = reshape(stimuliEnergy,[],res*res) * gaufun1(pp);

% Reshape back to images X orientations:
stimuliEnergyPRF = reshape(stimuliEnergyPRF,[],filt_prop.orientations);

%% 4. Calculate the standard deviation across orientation.

% Set a gain to scale the output:
gg = 1;

% Define a function that calculates the standard deviation across
% orientation:
modelfun = @(gg,dd) gg * std(dd,[],2);

modelOutput = modelfun(1,stimuliEnergyPRF);
