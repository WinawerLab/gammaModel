
% Filter many images with gabor patches and run through OV model and SOC
% model. Original code for images processing was developed by Kendrick Kay
% (PLOS Computational Biology, 2013).
%
% Images are from: - "Fred informed me that the citation for this database
% should be: Olmos, A., Kingdom, F. A. A. (2004), A biologically inspired
% algorithm for the recovery of shading and reflectance images, Perception,
% 33, 1463 - 1473. "
%
% Images are: 
% - converted to luminance values, convert to grayscale using [.3, .59, .11] for the RGB channels (rgb2gray.m) 
% - scaled .1 and 99.9 percentile to 0-1 
% - imresize 50% 
% - crop to square 
% - saved in the .mat file as "sqrt" units (uint8) 
% - 960 pixels x 960 pixels x 771 images.
%
% Dora Hermes, 2019

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
a1 = load(fullfile(dataDir,'stimuli','McGill_imageset','tempims.mat'));

% plot one image:
im0 = (double(a1.ims(:,:,50))/255).^2;  % now in luminance values between 0-1
 
% only do sqrt for viewing on your monitor (since your monitor performs a gamma squaring)
figure; imagesc(sqrt(im0));
colormap gray


%% %%%%%% preprocess images - part 1 is fast %%%%%%%%

%%%% DOWNSAMPLE TO INCREASE SPEED
% resize the stimuli to 240 x 240 to reduce computational time.
% use single-format to save memory.
temp = zeros(240,240,size(a1.ims,3),'single');
for pp = 1:size(a1.ims,3)
  temp(:,:,pp) = imresize(single(a1.ims(:,:,pp)),[240 240],'bilinear');
end
stimulus = temp;
clear temp;

%%%% RESCALE
% subtract off the background luminance (0.5).
% after these steps, all pixel values will be in the
% range [-.5,.5] with the background corresponding to 0.
stimulus = stimulus/255 - 0.5;

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
filt_prop.spacings = 1;             %   the separation of adjacent filters is 1 std dev of the Gaussian envelopes
                                    %     (this results in a 135 x 135 grid of positions)
filt_prop.orientations = 8;         %   filters occur at 8 orientations
filt_prop.phases = 2;               %   filters occur at 2 phases (between 0 and pi)
filt_prop.thres = 0.01;             %   the Gaussian envelopes are thresholded at .01
filt_prop.scaling = 2;              %   filters are scaled to have an equivalent Michelson contrast of 1
filt_prop.mode = 0;                 %   the dot-product between each filter and each stimulus is computed

% after this step, stimulus is images x phases*orientations*positions.
stimulus = applymultiscalegaborfilters(reshape(stimulus,270*270,[])', ...
	filt_prop.cycles,filt_prop.bandwidth,filt_prop.spacings,filt_prop.orientations,...
	filt_prop.phases,filt_prop.thres,filt_prop.scaling,filt_prop.mode);

save(fullfile(dataDir,'derivatives','gaborFilt','McGill_imageset_gaborFilt01.mat'),'stimulus','filt_prop')

%%
%% figure of images with different orientations

% load filtered textures used in paper
textures = load(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus');

% load natural images McGill set
load(fullfile(dataDir,'derivatives','gaborFilt','McGill_imageset_gaborFilt01.mat'),'stimulus','filt_prop')

% Calculate quadrature pairs: compute the square root of the sum of the
% squares of the outputs of quadrature-phase filter pairs (this is the
% standard complex-cell energy model). After this step, stimulus is images
% x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2)); % 
textures.stimulus = sqrt(blob(textures.stimulus.^2,2,2)); 
% --> stimulus is the input to the OV model

% Resolution of natural image and textures is the same:
res = sqrt(size(stimulus,2)/filt_prop.orientations);

%%%% Get images ready for SOC model
% compute the population term in the divisive-normalization equation.
% this term is simply the average across the complex-cell outputs
% at each position (averaging across orientation).
stimulusPOP = blob(stimulus,2,8)/8;
textures.stimulusPOP = blob(textures.stimulus,2,8)/8;

% Repeat the population term for each of the orientations
stimulusPOP = upsamplematrix(stimulusPOP,8,2,[],'nearest');
textures.stimulusPOP = upsamplematrix(textures.stimulusPOP,8,2,[],'nearest');

% Apply divisive normalization to the complex-cell outputs.  there are two
% parameters that influence this operation: an exponent term (r) and a
% semi-saturation term (s). the parameter values specified here were
% determined through a separate fitting procedure (see paper for details).
% for the purposes of this script, we will simply hard-code the parameter
% values here and not worry about attempting to fit the parameters.
r = 1;
s = 0.5;
stimulusSOC = stimulus.^r ./ (s.^r + stimulusPOP.^r);
textures.stimulusSOC = textures.stimulus.^r ./ (s.^r + textures.stimulusPOP.^r);
clear stimulusPOP;

% Sum across orientation. After this step, stimulusSOC is images x positions.
stimulusSOC = blob(stimulusSOC,2,8);
textures.stimulusSOC = blob(textures.stimulusSOC,2,8);

%%
%% Run images through SOC and OV models
%%


%% Get SOC & OV parameters all electrodes
%%%%% Subjects/electrodes:
subject_ind = [19 19 19 19 19 19  ... % S1
    24 24 ... % S2
    1001 1001 1001 1001 1001 1001 1001 1001]; % S3
electrodes = [107 108 109 115 120 121 ... % S1
    45 46 ... % S2
    49 50 52 57 58 59 60]; % S3
% define parameters to collect looping through electrodes/subjects:
socParams_all = zeros(length(electrodes),6);
ovParams_all = zeros(length(electrodes),5);
for ll = 1:length(electrodes)

    subj = subject_ind(ll);
    elec = electrodes(ll);
    
    % model fit on data
    analysisType = 'spectra200';
    
    % load SOC model fit
    modelType = 'fitSOCbbpower2';    
    load(fullfile(dataDir,'derivatives','gaborFilt','fitSOCbb',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'cross_SOCparams')
    % get median model parameters
    cross_SOCparams(cross_SOCparams(:,6)<0,6) = 0; % restrictrange at 0
    cross_SOCparams(cross_SOCparams(:,6)>1,6) = 1; % restrictrange at 1
    cross_SOCparams(:,3) = abs(cross_SOCparams(:,3)); % size>0
    socParams_all(ll,:) = median(cross_SOCparams);
    
    % load OV model fit
    modelType = 'OVsimple';   
    ov_exponents = [.1:.1:1];
    load(fullfile(dataDir,'derivatives','gaborFilt','deriveOV',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'cross_OVparams')   
    % get median model parameters and plot prediction
    ovParams_all(ll,:) = median(cross_OVparams(:,:,ov_exponents==.5));
end

%% Run all natural images through SOC and OV models
% Loop through stimuli and electrodes
% define output (stimuli X electrodes)
SOC_estimates = zeros(size(stimulus,1),size(socParams_all,1));
OV_estimates = zeros(size(stimulus,1),size(socParams_all,1));
for ss = 1:size(stimulus,1) % stimuli
    for ll = 1:size(socParams_all,1)
        [~,socEstimate] = helpfit_SOC2(stimulusSOC(ss,:),socParams_all(ll,:),[],[]);
        SOC_estimates(ss,ll) = socEstimate;
        [~,ovEstimate] = helpfit_OVexp(stimulus(ss,:),ovParams_all(ll,:),[],[]);
        OV_estimates(ss,ll) = ovEstimate;
    end
end

%% Run all textures images through SOC and OV models
% Loop through stimuli and electrodes
% define output (stimuli X electrodes)
textures.SOC_estimates = zeros(size(textures.stimulus,1),size(socParams_all,1));
textures.OV_estimates = zeros(size(textures.stimulus,1),size(socParams_all,1));
for ss = 1:size(textures.stimulus,1) % stimuli
    for ll = 1:size(socParams_all,1)
        [~,socEstimate] = helpfit_SOC2(textures.stimulusSOC(ss,:),socParams_all(ll,:),[],[]);
        textures.SOC_estimates(ss,ll) = socEstimate;
        [~,ovEstimate] = helpfit_OVexp(textures.stimulus(ss,:),ovParams_all(ll,:),[],[]);
        textures.OV_estimates(ss,ll) = ovEstimate;
    end
end

%%  Figure with output values directly from electrode estimates

figure('Position',[0 0 600 300])
subplot(1,2,1),hold on
plot(SOC_estimates,OV_estimates,'o','Color',[.5 .5 .5],'MarkerSize',2)

% add gratings:
plot(textures.SOC_estimates(39:46,:),textures.OV_estimates(39:46,:),'.','Color',[1 0 0])
plot(textures.SOC_estimates(50,:),textures.OV_estimates(50,:),'.','Color',[1 .5 .5])
xlim([min([textures.SOC_estimates(:);textures.OV_estimates(:)])-30 max([textures.SOC_estimates(:);textures.OV_estimates(:)])+30])
ylim([min([textures.SOC_estimates(:);textures.OV_estimates(:)])-30 max([textures.SOC_estimates(:);textures.OV_estimates(:)])+30])
axis square 

xlabel('bb prediction (SOC model)')
ylabel('gamma prediction (OV model)')
title('Natural Images(black) and Gratings(red)')

% set axis on limits from natural images:
subplot(1,2,2),hold on
plot(SOC_estimates,OV_estimates,'o','Color',[.5 .5 .5],'MarkerSize',2)

xlim([min([SOC_estimates(:);OV_estimates(:)])-10 max([SOC_estimates(:);OV_estimates(:)])+10])
ylim([min([SOC_estimates(:);OV_estimates(:)])-10 max([SOC_estimates(:);OV_estimates(:)])+10])
axis square 

xlabel('bb prediction (SOC model)')
ylabel('gamma prediction (OV model)')
title('Natural Images')

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',fullfile(dataDir,'derivatives','gaborFilt',...
%         ['NaturalImages_Test01']))
% print('-dpng','-r300',fullfile(dataDir,'derivatives','gaborFilt',...
%         ['NaturalImages_Test01']))

%% Figure normalized to grating response:

SOC_estimates_norm = zeros(size(SOC_estimates));
OV_estimates_norm = zeros(size(OV_estimates));
textures.SOC_estimates_norm = zeros(size(textures.SOC_estimates));
textures.OV_estimates_norm = zeros(size(textures.OV_estimates));

for kk = 1:size(SOC_estimates,2) % electrodes
    grat_SOC = textures.SOC_estimates(39,kk);
    SOC_estimates_norm(:,kk) = SOC_estimates(:,kk)./grat_SOC;
    textures.SOC_estimates_norm(:,kk) = textures.SOC_estimates(:,kk)./grat_SOC;
    
    grat_OV = textures.OV_estimates(39,kk);
    OV_estimates_norm(:,kk) = OV_estimates(:,kk)./grat_OV;
    textures.OV_estimates_norm(:,kk) = textures.OV_estimates(:,kk)./grat_OV;
end

figure('Position',[0 0 800 350])
subplot(4,4,[1:3 5:7 9:11])
hold on
plot(SOC_estimates_norm,OV_estimates_norm,'o','Color',[.5 .5 .5],'MarkerSize',2)
plot(textures.SOC_estimates_norm(39:46,:),textures.OV_estimates_norm(39:46,:),'o','Color',[1 .5 .5],'MarkerSize',2)
% [298,10] image/electrode example high OV
plot(SOC_estimates_norm(298,10),OV_estimates_norm(298,10),'o','Color',[.2 .2 1],'MarkerSize',2)
% [655,1] image/electrode example high SOC
plot(SOC_estimates_norm(655,1),OV_estimates_norm(655,1),'o','Color',[.2 1 .2],'MarkerSize',2)

xlim([min([SOC_estimates_norm(:);OV_estimates_norm(:)]) max([SOC_estimates_norm(:);OV_estimates_norm(:)])])
ylim([min([SOC_estimates_norm(:);OV_estimates_norm(:)]) max([OV_estimates_norm(:)]+1)])
% ylim([min([SOC_estimates_norm(:);OV_estimates_norm(:)]) 8])
plot([0 1],[1 1],'k','LineWidth',1)
plot([1 1],[0 1],'k','LineWidth',1)

ylabel('gamma prediction (OV model)')
title('Natural images normalized to grating(red)')
set(gca,'XTick',[1 10 20],'YTick',[1])

subplot(4,4,[4 8 12])
[nn,xx] = hist(OV_estimates_norm(:),[0:.05:1.1]);
% barh(xx,nn,'FaceColor',[.5 .5 .5]);
plot(nn,xx,'Color',[.5 .5 .5]);
ylim([min([SOC_estimates_norm(:);OV_estimates_norm(:)]) max([OV_estimates_norm(:)]+1)])
xlim([0 max(nn)+10])
hold on
% all_grating_OV = textures.OV_estimates_norm(39:46,:);
% [nn,xx] = hist(all_grating_OV(:),[0:.05:1.1]);
% % barh(xx,nn,'FaceColor',[1 .5 .5]);
% plot(nn,xx,'Color',[1 .5 .5]);
set(gca,'YTick',[1])
box off

subplot(4,4,13:15)
[nn,xx] = hist(SOC_estimates_norm(:),[0:.05:20]);
% bar(xx,nn,'FaceColor',[.5 .5 .5]);
plot(xx,nn,'Color',[.5 .5 .5],'LineWidth',1);
xlim([min([SOC_estimates_norm(:);OV_estimates_norm(:)]) max([SOC_estimates_norm(:);OV_estimates_norm(:)])])
ylim([0 max(nn)+10])
hold on
% all_grating_SOC = textures.SOC_estimates_norm(39:46,:);
% [nn,xx] = hist(all_grating_SOC(:),[0:.05:20]);
% bar(xx,nn,'FaceColor',[1 .5 .5]);
% plot(xx,nn,'Color',[1 .5 .5],'LineWidth',1);
box off
set(gca,'XTick',[1 10 20])
xlabel('bb prediction (SOC model)')

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',fullfile(dataDir,'derivatives','gaborFilt',...
%         ['NaturalImages_Test03_norm']))
% print('-dpng','-r300',fullfile(dataDir,'derivatives','gaborFilt',...
%         ['NaturalImages_Test03_norm']))

%% Display some images with large normalized values:
%%
[ii,jj] = find(OV_estimates_norm>.200);
ii = ii(24);
jj= jj(24);
% [298,10] image/electrode

% plot one image:
im0 = (double(a1.ims(:,:,ii))/255).^2;  % now in luminance values between 0-1 
% only do sqrt for viewing on your monitor (since your monitor performs a gamma squaring)
figure; imagesc(sqrt(im0));
colormap gray

seed_params = ovParams_all(jj,:);
figure('Position',[0 0 500 500])
imagesc(sqrt(im0)),hold on
% Move pRF in 135x135 to original size of 960x960 
%   half of the screen is 60 in filtered and 960 in original
orig_x = 480 + (seed_params(1)-res/2) * 480./60;
orig_y = 480 + (seed_params(2)-res/2) * 480./60;
orig_sigma = seed_params(3) * 480./60;
axis image
colormap gray

numPoints = 50;
c.th = linspace(0,2*pi, numPoints);
[c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*orig_sigma);
plot(c.x + orig_y, c.y + orig_x, 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
[c.x, c.y] = pol2cart(c.th, 2*ones(1,numPoints)*orig_sigma);
plot(c.x + orig_y, c.y + orig_x, 'r:') % this is just reversed because plot and imagesc are opposite, checked this with contour

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'derivatives','gaborFilt',...
        ['NaturalImages_Test02_example']))
print('-dpng','-r300',fullfile(dataDir,'derivatives','gaborFilt',...
        ['NaturalImages_Test02_example']))

%% Display some images with large normalized SOC values:
%%    
[ii,jj] = find(SOC_estimates_norm>20);
ii = ii(3);
jj= jj(3);
% [655,1] image/electrode

% plot one image:
im0 = (double(a1.ims(:,:,ii))/255).^2;  % now in luminance values between 0-1 
% only do sqrt for viewing on your monitor (since your monitor performs a gamma squaring)
figure; imagesc(sqrt(im0));
colormap gray

seed_params = ovParams_all(jj,:);
figure('Position',[0 0 500 500])
imagesc(sqrt(im0)),hold on
% Move pRF in 135x135 to original size of 960x960 
%   half of the screen is 60 in filtered and 960 in original
orig_x = 480 + (seed_params(1)-res/2) * 480./60;
orig_y = 480 + (seed_params(2)-res/2) * 480./60;
orig_sigma = seed_params(3) * 480./60;
axis image
colormap gray

numPoints = 50;
c.th = linspace(0,2*pi, numPoints);
[c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*orig_sigma);
plot(c.x + orig_y, c.y + orig_x, 'r') % this is just reversed because plot and imagesc are opposite, checked this with contour
[c.x, c.y] = pol2cart(c.th, 2*ones(1,numPoints)*orig_sigma);
plot(c.x + orig_y, c.y + orig_x, 'r:') % this is just reversed because plot and imagesc are opposite, checked this with contour

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'derivatives','gaborFilt',...
        ['NaturalImages_Test02_example2']))
print('-dpng','-r300',fullfile(dataDir,'derivatives','gaborFilt',...
        ['NaturalImages_Test02_example2']))

%% Display some images with large values:
%%
[ii,jj] = find(OV_estimates>600);

figure
imagesc(reshape(stimulusSOC(ii,:),res,res))

% plot one image:
im0 = (double(a1.ims(:,:,ii))/255).^2;  % now in luminance values between 0-1 
% only do sqrt for viewing on your monitor (since your monitor performs a gamma squaring)
figure; imagesc(sqrt(im0));
colormap gray





%%

ims_plot = [689];%[145 146 342 493 689];

figure
for kk = 1:length(ims_plot)
    % plot one image:
    im0 = (double(a1.ims(:,:,ims_plot(kk)))/255).^2;  % now in luminance values between 0-1 
    % only do sqrt for viewing on your monitor (since your monitor performs a gamma squaring)
    subplot(1,length(ims_plot),kk)
    imagesc(sqrt(im0));
    colormap gray
    axis image
    axis off
end

%%
%% Figure with output from 1 electrode
%%

%%  Figure with output values directly from electrode estimates

example_elec1 = 15;
% example_elec2 = 10;

figure('Position',[0 0 1200 800])

for example_elec1 = 1:15
    subplot(3,5,example_elec1),hold on
    h = scatter(SOC_estimates(:,example_elec1),OV_estimates(:,example_elec1),15,[.2 .2 .2],'filled','MarkerFaceAlpha',.3);
    % add gratings:
    scatter(textures.SOC_estimates(39:46,example_elec1),textures.OV_estimates(39:46,example_elec1),15,[.5 0 0],'filled')
    scatter(textures.SOC_estimates(50,example_elec1),textures.OV_estimates(50,example_elec1),15,[.8 0 0],'filled')
    scatter(textures.SOC_estimates(49,example_elec1),textures.OV_estimates(49,example_elec1),15,[1 0 0],'filled')
    scatter(textures.SOC_estimates(48,example_elec1),textures.OV_estimates(48,example_elec1),15,[1 .2 .2],'filled')
    scatter(textures.SOC_estimates(47,example_elec1),textures.OV_estimates(47,example_elec1),15,[1 .4 .4],'filled')
    axis square 
    axis tight

    xlabel('SOC prediction')
    ylabel('OV prediction')
    title(['El' int2str(example_elec1) ' Natural (gray) Gratings (red)'])
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'derivatives','gaborFilt',...
        ['NaturalImages_Test03']))
print('-dpng','-r300',fullfile(dataDir,'derivatives','gaborFilt',...
        ['NaturalImages_Test03']))

