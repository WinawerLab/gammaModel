
%
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
load(fullfile(dataDir,'soc_bids','stimuli','task-soc_stimuli.mat'),'stimuli')


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

save(fullfile(dataDir,'soc_bids','derivatives','stimuliGaborfilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

%% %%%%%% preprocess images - part 3 is fast %%%%%%%%

load(fullfile(dataDir,'soc_bids','derivatives','stimuliGaborfilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

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
save(fullfile(dataDir,'soc_bids','derivatives','stimuliGaborfilt','task-soc_stimuli_gaborFilt02.mat'),'stimulus')

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
load(fullfile(dataDir,'soc_bids','derivatives','stimuliGaborfilt','task-soc_stimuli_gaborFilt02.mat'),'stimulus')


%%
%% load ECoG data 
%%

dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
%%%%% Pick a subject:
subjects = [19,23,24];
s = 1; subj = subjects(s);

%%%%% Pick an electrode:
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3
elec = 109;

% Choose an analysis type:
% analysisType = 'spectraRERP500';
% analysisType = 'spectra';
analysisType = 'spectra200';

% load the data:
dataFitName = [dataRootPath '/sub-' int2str(subj) '/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
load(dataFitName,'resamp_parms')

% Broadband estimate, one value per image:
bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
ecog_bb = squeeze(median(resamp_parms(:,:,2),2))-bb_base;
ecog_bb_err=[squeeze(quantile(resamp_parms(:,:,2),.025,2)) ...
    squeeze(quantile(resamp_parms(:,:,2),.975,2))]'-bb_base;

ecog_g = squeeze(median(resamp_parms(:,:,3),2));
ecog_g_err=[squeeze(quantile(resamp_parms(:,:,3),.025,2)) ...
    squeeze(quantile(resamp_parms(:,:,3),.975,2))]';

ecog_a = squeeze(median(resamp_parms(:,:,6),2));
ecog_a_err=[squeeze(quantile(resamp_parms(:,:,6),.025,2)) ...
    squeeze(quantile(resamp_parms(:,:,6),.975,2))]';

% Pick ECoG values to analyse:
ecog_vals=ecog_bb;

%% %%% Prepare for model fitting

% to prepare for the call to fitnonlinearmodel.m, we have to define various
% input parameters.  this is what we will now do.

% define constants
res = 135;  % resolution of the pre-processed stimuli

% in this example script, we will be fitting only some of the parameters of the
% SOC model.  the other parameters have already been fixed (see above).
%
% the parameters that we will be fitting are [R C S G N C] where
%   R is the row index of the center of the 2D Gaussian
%   C is the column index of the center of the 2D Gaussian
%   S is the standard deviation of the 2D Gaussian
%   G is a gain parameter
%   N is the exponent of the power-law nonlinearity
%   C is a parameter that controls the strength of second-order contrast

% define initial seeds for model parameters.  to help avoid the problem of
% local minima, we will start at several different initial seeds (corresponding to
% different combinations of N and C), and we will choose the model that
% achieves the best fit to the data.
Ns = [.05 .1 .2 .3 .4 .5 .6 .7 1];
Cs = [.1 .4 .7 .8 .85 .9 .95 .975 .99 .995];
seeds = [];
for p=1:length(Ns)
  for q=1:length(Cs)
    seeds = cat(1,seeds,[(1+res)/2 (1+res)/2 res/4*sqrt(0.5) 10 Ns(p) Cs(q)]);  
  end
end

% define bounds for the model parameters
bounds = [1-res+1 1-res+1 0   -Inf 0   0;
          2*res-1 2*res-1 Inf  Inf Inf 1];

% define a version of bounds where we insert NaNs in the first row
% in the spots corresponding to the N and C parameters.  this
% indicates to fix these parameters and not optimize them.
boundsFIX = bounds;
boundsFIX(1,5:6) = NaN;

% issue a dummy call to makegaussian2d.m to pre-compute xx and yy.
% these variables are re-used to achieve faster computation.
[d,xx,yy] = makegaussian2d(res,2,2,2,2);

% define a helper function that we will use for the SOC model.
% this function accepts stimuli (dd, a matrix of size A x 90*90),
% weights (wts, a vector of size 135*135 x 1), and a parameter
% (c, a scalar) and outputs a measure of variance (as a vector
% of size A x 1).  intuitively, this function computes
% a weighted average, subtracts that off, squares, and
% computes a weighted sum.
socfun = @(dd,wts,c) bsxfun(@minus,dd,c*(dd*wts)).^2 * wts;

% define another helper function.  given a set of parameters,
% this function outputs a 2D Gaussian (as a vector of size 135*135 x 1)
gaufun = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

% we will now define a function that implements the SOC model.  this function
% accepts a set of parameters (pp, a vector of size 1 x 6) and a set of stimuli
% (dd, a matrix of size A x 135*135) and outputs the predicted response to those
% stimuli (as a vector of size A x 1).
modelfun = @(pp,dd) pp(4)*(socfun(dd,gaufun(pp),restrictrange(pp(6),0,1)).^pp(5));

% notice that in the above function, we use restrictrange to force C to be
% between 0 and 1.  this is necessary because even though we defined bounds
% for the parameters (see above), we will be using the Levenberg-Marquardt
% optimization algorithm, which does not respect parameter bounds.

% we are ready to define the final model specification.  in the following,
% we specify a stepwise fitting scheme.  in the first fit (the first row),
% we start at the seed and optimize all parameters except the N and C parameters.
% in the second fit (the second row), we start at the parameters estimated in
% the first fit and optimize all parameters.  the reason that the first entry
% is [] is that we will be using a mechanism that evaluates multiple initial
% seeds, and in that case the first entry here is ignored.
model = {{[]         boundsFIX   modelfun} ...
         {@(ss) ss   bounds      @(ss) modelfun}};

% define optimization options (these are added on top of
% the default options used in fitnonlinearmodel.m)
optimoptions = {'Algorithm' 'levenberg-marquardt' 'Display' 'off'};

% define the resampling scheme to use.  here, we use 0, which
% means to just fit the data (no cross-validation nor bootstrapping).
resampling = 0;

% define the metric that we will use to quantify goodness-of-fit.
metric = @(a,b) calccod(a,b,[],[],0);

% specify the index of the voxel to fit
% ix = 204; 

% finally, construct the options struct that will be passed to fitnonlinearmodel.m
opt = struct( ...
  'stimulus',     stimulus, ...
  'data',         ecog_vals, ...
  'model',        {model}, ...
  'seed',         seeds, ...
  'optimoptions', {optimoptions}, ...
  'resampling',   resampling, ...
  'metric',       metric);

% do a quick inspection of opt
opt


%% %%% Fit the model

results = fitnonlinearmodel(opt);


%% %%% Inspect the results

% these are the final parameter estimates
results.params
% the fitted parameters are [R C S G N C] where
%   R is the row index of the center of the 2D Gaussian
%   C is the column index of the center of the 2D Gaussian
%   S is the standard deviation of the 2D Gaussian
%   G is a gain parameter
%   N is the exponent of the power-law nonlinearity
%   C is a parameter that controls the strength of second-order contrast

% this is the R^2 between the model fit and the data
results.trainperformance
%%
% visualize the data and the model fit
figure; setfigurepos([100 100 600 450]); 

subplot(2,1,1),hold on;
bar(ecog_vals,1);
% errorbar2(1:86,betamn(ix,1:99),betase(ix,1:99),'v','g-','LineWidth',2);
modelfit = modelfun(results.params,stimulus);

plot(1:86,modelfit,'r-','LineWidth',3);

% plot stimulus cutoffs
stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
for k=1:length(stim_change)
    plot([stim_change(k) stim_change(k)],[0 1],'g')
end
text([19 40 46 52 55 61 70 74 79 84],zeros(10,1)-.25,{'space','orie','grat','pl','circ','zcon','sp','zsp','coh','nm'})

ax = axis;
axis([0 100 ax(3:4)]);
xlabel('Stimulus number');
ylabel('ECoG signal');
title(['Data and model fit R=' int2str(results.trainperformance)]);

subplot(2,1,2)
%%% LOOK AT WHERE THE GAUSSIAN IS
gau = makegaussian2d(res,results.params(1),results.params(2),results.params(3)/sqrt(results.params(5)),results.params(3)/sqrt(results.params(5)),xx,yy,0,0);
imagesc(gau);axis equal; colorbar;hold on

% set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/fit_tests/' subj '_el'  int2str(electr) '_gainseed_' int2str(seeds(1,4))])
% print('-depsc','-r300',['./figures/fit_tests/' subj '_el'  int2str(electr) '_gainseed_' int2str(seeds(1,4))])

