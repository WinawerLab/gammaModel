%
%
%
% This scripts is used to refit fMRI data with the SOC model
%
% !!! Not sure how the images should be filtered and fitted here, should
% redo for paper
%
%

clear all

% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

%% Load fMRI data
%
% load fMRI stimulus numbers that were used in ECoG
load(fullfile(dataDir,'bold','imnumbersbold2ecog')); % loads imnumbers_ecog

% load fMRI beta values
load(fullfile(dataDir,'bold','dataset03.mat'),'betamn','betase','roi','vxsselect','roilabels');
% roi: roi label numbers for each voxel
% roilabels: names corresponding to roi numbers

% load modelfit to fMRI beta values
a = load(fullfile(dataDir,'bold','dataset03_benchmark.mat'));

% load all fMRI stimuli
load(fullfile(dataDir,'bold','stimuli.mat'))


%% xys is not saved with Kay paper and we need to refit the SOC model 

% Note that these are filtered stimuli, not sure how this was done...
load('./data/stimulus_filt','stimulus')

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

% each of the 99 stimuli that we are considering actually consists of 9 frames.
% we will be handling this by computing, for each stimulus, the predicted response
% for each of the 9 frames and then averaging across these predicted responses.
% the legwork for doing this is handled by the fitting function we will be using,
% fitnonlinearmodel.m.  that function requires that the individual frames that comprise
% a single stimulus be arranged along the third dimension, so that is now what we do.

% reshape stimuli into desired format: 99 stimuli x 90*90 positions x 9 frames.
stimulus = permute(reshape(stimulus,9,77,8100),[2 3 1]);

% inspect one of the stimuli
figure;
mx = max(abs(stimulus(:)));
imagesc(reshape(stimulus(12,:,1),[90 90]));
axis image tight;
caxis([0 mx]);
colormap(gray);
colorbar;
title('Stimulus');

%%

% loop over all V1/V2/V3-voxels
voxels2fit = find(ismember(roi,[1 2 3 4]) & vxsselect(:,3)==1 & benchmark_perf.modelpredr2>.5);

for v_ind=1:length(voxels2fit)

    ix=voxels2fit(v_ind);

    % to prepare for the call to fitnonlinearmodel.m, we have to define various
    % input parameters.  this is what we will now do.

    % define constants
    res = 90;  % resolution of the pre-processed stimuli

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
    Ns = [.05 .1 .3 .5];
    Cs = [.4 .7 .9 .95];
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
    % weights (wts, a vector of size 90*90 x 1), and a parameter
    % (c, a scalar) and outputs a measure of variance (as a vector
    % of size A x 1).  intuitively, this function computes
    % a weighted average, subtracts that off, squares, and
    % computes a weighted sum.
    socfun = @(dd,wts,c) bsxfun(@minus,dd,c*(dd*wts)).^2 * wts;

    % define another helper function.  given a set of parameters,
    % this function outputs a 2D Gaussian (as a vector of size 90*90 x 1)
    gaufun = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

    % we will now define a function that implements the SOC model.  this function
    % accepts a set of parameters (pp, a vector of size 1 x 6) and a set of stimuli
    % (dd, a matrix of size A x 90*90) and outputs the predicted response to those
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
    % metric = @(a,b) calccod(a,b,[],[],0);
    metric = @(a,b) calccod(a,b);

    % finally, construct the options struct that will be passed to fitnonlinearmodel.m
    opt = struct( ...
      'stimulus',     stimulus(1:38,:,:), ...
      'data',         betamn(ix,imnumbers_ecog(1:38))', ...
      'model',        {model}, ...
      'seed',         seeds, ...
      'optimoptions', {optimoptions}, ...
      'resampling',   resampling, ...
      'metric',       metric);

    % do a quick inspection of opt
    % opt

    results = fitnonlinearmodel(opt);

    save(['./data/fit_spatial_bold/voxel_' int2str(ix)],'results','ix');

    disp(['done voxel ' int2str(v_ind) ' out of ' length(voxels2fit)])

end

save(['./data/fit_spatial_bold/opt_fitspatial'],'opt');

%% Save fit for voxels with good fmri results in V1-V4

clear all

% Load data from the third dataset
load('../../bold/dataset03.mat','betamn','betase','roi','vxsselect','roilabels');
% Load model performance, to only get voxels with pretty good fit
benchmark_perf = load('../../bold/dataset03_benchmark.mat');
% Load image numbers to use
load('./data/imnumbersbold2ecog','imnumbers_ecog')

voxels2fit = find(ismember(roi,[1 2 3 4]) & vxsselect(:,3)==1 & benchmark_perf.modelpredr2>.5);

%%%%% get prf in degrees for fMRI
x_y_fmri=NaN(length(roi),2);

for v_ind=1:length(voxels2fit)

    ix=voxels2fit(v_ind);

    load(['./data/fit_spatial_bold/voxel_' int2str(ix)],'results','ix');
    load(['./data/fit_spatial_bold/opt_fitspatial'],'opt');

    % convert prf estimate to an estimate degrees of visual angle
    % BOLD: 
    %   images  150 pixels, 12.5 degrees
    %   resized 180 pixels, 15   degrees, 90  filter outputs
    % ECoG: 
    %   images  240 pixels, 20   degrees
    %   resized 270 pixels, 22.5 degrees, 135 filter outputs

    x_y_fmri(ix,:)=(results.params(1:2)-45)*15/90; % subtract center and convert to degrees
end

save('./data/fit_prf_location_BOLD','x_y_fmri')



%% Match each electrode to some voxels:
% clear all

load('./data/fit_prf_location_BOLD','x_y_fmri')

load('../../bold/dataset03.mat','betamn','betase','roi','vxsselect','roilabels');
load('./data/imnumbersbold2ecog','imnumbers_ecog')
clear out
subj_nr = 19;
electrodes = [108 109 115 120 121 107];

for el = 5%1:length(electrodes)
    electr = electrodes(el);

    [v_area,xys,roi_labels] = subj_prf_info(subj_nr,electr);
    
    out(el).el = electr;
    out(el).area = v_area;
    out(el).xys = xys; % ECoG pRF experiment x, y, size in degrees of visual angle 
    
    % find voxels closest to this pRF x and y 
    prf_dist_fmri_ecog = sqrt(sum((x_y_fmri-repmat(xys(1:2),size(x_y_fmri,1),1)).^2,2));
    
    % exclude voxels from non-roi areas:
    prf_dist_fmri_ecog(~ismember(roi,v_area)) = NaN;

    [vox_ind,vox_near] = sort(prf_dist_fmri_ecog);

    % Take the closest 5 voxels near an ROI, this distance is in degree of
    % visual angle
    out(el).voxels_match = vox_near(1:5);
    out(el).distance = vox_ind(1:5);
    
    % TODO: remove voxels if the difference in the degree of visual angle is too large
    
    out(el).voxels_roi = roi(out(el).voxels_match);
    out(el).betamatch = betamn(out(el).voxels_match,imnumbers_ecog);
    
    % save out
%     save(['./data/match_ECoG_BOLD_subj' int2str(subj_nr) '_elec' int2str(electr)],'out')
end
%%
figure
for el = 1:length(electrodes)
    subplot(length(electrodes),1,el),hold on
    plot(mean(out(el).betamatch))
    title(int2str(electrodes(el)))
    axis tight
end  
