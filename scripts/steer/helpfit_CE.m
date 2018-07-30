function [results,modelfit] = helpfit_CE(stimulus,params,ecog_vals,seeds)
%
% helpfunction to fit or apply the CE (contrast energy) model
%
% inputs
%       stimulus: the contrast images
%       params: parameters for the CE model 
%               if ~isempty - apply model without fitting 
%               if empty - fit model 
%             parameters are [R C S G] where:
%             R is the row index of the center of the 2D Gaussian
%             C is the column index of the center of the 2D Gaussian
%             S is the standard deviation of the 2D Gaussian
%             G is a gain parameter
%       ecog_vals: values to fit
%       seeds: seeds for fitting, leave empty if only applying model 
%
% output
%       results: 
%       modelfit: model applied to stimulus
%
% example:
%       [results,modelfit] = helpfit_CE(stimulus,params,ecog_vals,seeds)
%

results = [];

% get resolution of the stimulus
res = sqrt(size(stimulus,2));

% issue a dummy call to makegaussian2d.m to pre-compute xx and yy.
% these variables are re-used to achieve faster computation.
[~,xx,yy] = makegaussian2d(res,2,2,2,2);

% define another helper function.  given a set of parameters,
% this function outputs a 2D Gaussian (as a vector of size 135*135 x 1)
gaufun = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

% we will now define a function that implements the SOC model.  this function
% accepts a set of parameters (pp, a vector of size 1 x 6) and a set of stimuli
% (dd, a matrix of size A x 135*135) and outputs the predicted response to those
% stimuli (as a vector of size A x 1).
modelfun = @(pp,dd) pp(4)*dd*gaufun(pp);

% notice that in the above function, we use restrictrange to force C to be
% between 0 and 1.  this is necessary because even though we defined the
% following bounds for the parameters (see below), we will be using the
% Levenberg-Marquardt optimization algorithm, which does not respect
% parameter bounds. Define bounds for the model parameters:
bounds = [1-res+1 1-res+1 1   -Inf;
          2*res-1 2*res-1 Inf  Inf];

% we are ready to define the final model specification.  in the following,
% we specify a stepwise fitting scheme.  in the first fit (the first row),
% we start at the seed and optimize all parameters except the N and C parameters.
% in the second fit (the second row), we start at the parameters estimated in
% the first fit and optimize all parameters.  the reason that the first entry
% is [] is that we will be using a mechanism that evaluates multiple initial
% seeds, and in that case the first entry here is ignored.
model = {{[]         bounds   modelfun}};

% define optimization options (these are added on top of
% the default options used in fitnonlinearmodel.m)
optimoptions = {'Algorithm' 'levenberg-marquardt' 'Display' 'iter'};

% define the resampling scheme to use.  here, we use 0, which
% means to just fit the data (no cross-validation nor bootstrapping).
resampling = 0;

% define the metric that we will use to quantify goodness-of-fit.
metric = @(a,b) calccod(a,b);
% metric = @(a,b) calccod(a,b,[],[],0);

if isempty(params) % do the fitting
    % set the fitting inputs:
    opt = struct( ...
      'stimulus',     stimulus, ...
      'data',         ecog_vals, ...
      'model',        {model}, ...
      'seed',         seeds, ...
      'optimoptions', {optimoptions}, ...
      'resampling',   resampling, ...
      'metric',       metric);

    % do the fitting:
    results = fitnonlinearmodel(opt);
    % return the estimates:
    modelfit = modelfun(results.params,stimulus);
    
else % if parameters are there just apply the model with these params
    % return the estimates:
    modelfit = modelfun(params,stimulus);
end
