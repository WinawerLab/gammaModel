
% Filter patterns with gabor patches. Original code was developed by
% Kendrick Kay (PLOS Computational Biology), Code is used here to filter
% patterns from ECoG experiments.
%
% Dora Hermes, 2017

clear all
% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath('~/Documents/git/ecogBasicCode/')
addpath(genpath('~/Documents/m-files/knkutils'));

%%
%% %%%%%% START preprocess images %%%%%%%%
%%

% load images
load(fullfile(dataDir,'stimuli','task-soc_stimuli.mat'),'stimuli')


%% %%%%%% preprocess images - part 1 is fast %%%%%%%%

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
% pad the stimulus with zeros (to reduce edge effects).
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

%% %%%%%% preprocess images - part 2 is timeconsuming %%%%%%%%

% Apply Gabor filters to the stimuli.  filters occur at different positions,
% orientations, and phases.  there are several parameters that govern the
% design of the filters:
filt_prop.cycles = 60*(270/240);    %   the number of cycles per image is 60*(270/240)
filt_prop.bandwidth = -1;           %   the spatial frequency bandwidth of the filters is 1 octave
filt_prop.spacings=1;               %   the separation of adjacent filters is 1 std dev of the Gaussian envelopes
                                    %     (this results in a 90 x 90 grid of positions)
filt_prop.orientations=8;           %   filters occur at 8 orientations
filt_prop.phases=2;                 %   filters occur at 2 phases (between 0 and pi)
filt_prop.thres=0.01;               %   the Gaussian envelopes are thresholded at .01
filt_prop.scaling=2;                %   filters are scaled to have an equivalent Michelson contrast of 1
filt_prop.mode=0;                   %   the dot-product between each filter and each stimulus is computed

% after this step, stimulus is images x phases*orientations*positions.
stimulus = applymultiscalegaborfilters(reshape(stimulus,270*270,[])', ...
  filt_prop.cycles,filt_prop.bandwidth,filt_prop.spacings,filt_prop.orientations,...
  filt_prop.phases,filt_prop.thres,filt_prop.scaling,filt_prop.mode);

save(fullfile(dataDir,'derivatives','Gaborfilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

%% figure of images with different orientations

load(fullfile(dataDir,'derivatives','Gaborfilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')
stimulus = sqrt(blob(stimulus.^2,2,2)); % quadrature pairs

numOrientations = 8;

res = sqrt(size(stimulus,2)/numOrientations);


% Now to view an image, one can use:
inNr = 12;
thisIm = reshape(stimulus(inNr,:),numOrientations,res,res); 

% Make the figure
figure('Position',[0 0 1200 300]), 

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

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','filteredStimFigure',...
    ['stimulus-' int2str(inNr) '_filtered_orientations']))
print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','filteredStimFigure',...
    ['stimulus-' int2str(inNr) '_filtered_orientations']))

%% %%%%%% preprocess images - part 3 is fast %%%%%%%%

load(fullfile(dataDir,'soc_bids','derivatives','Gaborfilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

% compute the square root of the sum of the squares of the outputs of
% quadrature-phase filter pairs (this is the standard complex-cell energy model).
% after this step, stimulus is images x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2));

% compute the population term in the divisive-normalization equation.
% this term is simply the average across the complex-cell outputs
% at each position (averaging across orientation).
stimulusPOP = blob(stimulus,2,8)/8;

% repeat the population term for each of the orientations
stimulusPOP = upsamplematrix(stimulusPOP,8,2,[],'nearest');

% apply divisive normalization to the complex-cell outputs.  there are two parameters
% that influence this operation: an exponent term (r) and a semi-saturation term (s).
% the parameter values specified here were determined through a separate fitting
% procedure (see paper for details).  for the purposes of this script, we will
% simply hard-code the parameter values here and not worry about attempting to fit
% the parameters.
r = 1;
s = 0.5;
stimulus = stimulus.^r ./ (s.^r + stimulusPOP.^r);
clear stimulusPOP;

% sum across orientation.  after this step, stimulus is images x positions.
stimulus = blob(stimulus,2,8);
save(fullfile(dataDir,'derivatives','Gaborfilt','task-soc_stimuli_gaborFilt02.mat'),'stimulus')

%%
% inspect one of the stimuli
figure;
mx = max(abs(stimulus(:)));
imagesc(reshape(stimulus(12,:),[135 135]));
axis image tight;
caxis([0 mx]);
colormap(gray);
colorbar;
title('Stimulus');
    
%%
%% %%%%%% END preprocess images %%%%%%%%
%%


%%
%% %%%%%% LOAD preprocessed images  - skip timeconsuming part %%%%%%%%
%%
load(fullfile(dataDir,'derivatives','Gaborfilt','task-soc_stimuli_gaborFilt02.mat'),'stimulus')


%%
%% Visualize energy per orientation
%%

load(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')
stimulus = sqrt(blob(stimulus.^2,2,2)); % quadrature pairs
% sum across orientation.  after this step, stimulus is images x positions.
imEnergyMean = blob(stimulus,2,8);

numOrientations = 8;

res = sqrt(size(stimulus,2)/numOrientations);

figure('Position',[0 0 250 500]),
im_nrs = [50 49 54 10];

for ii = 1:length(im_nrs)
    % Select an image number:
    inNr = im_nrs(ii);
    thisIm = reshape(stimulus(inNr,:),numOrientations,res,res); 

    % Imagine this pRF
    pp = [res/2 res/2 res/5];
    [~,xx,yy] = makegaussian2d(res,2,2,2,2);
    gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
        pp(3),xx,yy,0,0)/(2*pi*pp(3)^2)); % Gaussian or pRF
    imEnergyPrf = imEnergyMean(inNr,:)'.*gaufun1(pp);

    % Make figure for energy for every orientation
    OrientationEnergy = zeros(1,numOrientations);
    for kk = 1:numOrientations
        thisImPrf = thisIm(kk,:,:);
        thisImPrf = thisImPrf(:).*gaufun1(pp);
        % Get energy for this orientation
        OrientationEnergy(kk) = sum(thisImPrf);
    end
    
    orient_order = [3 4 5 6 7 8 1 2];

    % Plot original image and energy across orientations
    subplot(length(im_nrs),2,ii*2-1)
    mx = max(abs(stimuli(:)));
    imagesc(stimuli(:,:,inNr))
    axis image tight off;
    caxis([0 mx]);
    colormap(gray);

    subplot(length(im_nrs),2,ii*2)
    bar(OrientationEnergy(orient_order),'w')
    xlim([0 9]),ylim([0 .2])
    box off
    ylabel('energy')
    title(['im ' int2str(inNr) ' sd=' num2str(var(OrientationEnergy).^.5,3)])
    set(gca,'FontName','Arial','FontSize',8)

end
 
% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',fullfile(dataDir,'derivatives','gaborFilt','filteredStimFigure',...
%     ['4stimuli_orientationEnergy']))
% print('-depsc','-r300',fullfile(dataDir,'derivatives','gaborFilt','filteredStimFigure',...
%     ['4stimuli_orientationEnergy']))





