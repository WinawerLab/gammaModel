
% Correlate image contrast of patterns filtered with steerable pyramid
% filters with gamma/broadband power.
%
% Dora Hermes, 2017

clear all
% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath('/Users/dora/Documents/git/ecogBasicCode/')
addpath(genpath('/Users/dora/Documents/m-files/knkutils'));
addpath('/Users/dora/Documents/m-files/steer/')

%% Load downsampled energies:

% load images
load(fullfile(dataDir,'soc_bids','stimuli','task-soc_stimuli.mat'),'stimuli')

% load downsamples image energies 
load(fullfile(dataDir,'soc_bids','derivatives','stimuliSteerfilt',['task-soc_stimuli_steerFilt01.mat']),'imEnergy')
% These filters have 5 levels:
% Level 1: 16 cycles per image: 1 cycle ~ 1.25 degrees
% Level 2: 8 cycles per image: 1 cycle ~ 2.5 degrees
% Level 3: 4 cycles per image: 1 cycle ~ 5 degrees
% Level 4: 2 cycles per image: 1 cycle ~ 10 degrees
% Level 5: 1 cycles per image: 1 cycle ~ 20 degrees

imEnergy_orig = imEnergy;
clear imEnergy;

%% Load ECoG data and fit

%%%%% Pick a subject:
subjects = [19,23,24];
s = 1; subj = subjects(s);

nrOrientations = 8;

if subj == 9
    im_deg = rad2deg(atan(20.7./61));
elseif subj == 19
    im_deg = rad2deg(atan(17.9./50));
    electrodes = [107 108 109 115 120 121]; % S1
elseif subj==23
    im_deg = rad2deg(atan(17.9./36));
    % electrodes = [53 54]; % S2
elseif subj==24
    im_deg = rad2deg(atan(17.9./45));
    electrodes = [45 46]; % S3
end
electrodes = [109];

res = sqrt(size(imEnergy_orig,2)/nrOrientations);  % resolution of the pre-processed stimuli

for el = 1:length(electrodes)
    elec = electrodes(el);

    [v_area,xys,roi_labels] = subj_prf_info(subj,elec);
    % Convert xys from degrees to pixels
    xys_pix = [res./2 res./2 0] + res.*(xys./im_deg);
    xys_pix(1:2) = [res-xys_pix(2) xys_pix(1)]; % make sure that it matches images

    if xys(3)<0.50 % small pRF size
        imEnergy = sum(imEnergy_orig(:,:,1),3);
    elseif xys(3)>=0.50 && xys(3)< 0.9 
       imEnergy = sum(imEnergy_orig(:,:,1),3);
    elseif xys(3)>=0.9 && xys(3)< 1.7 
       imEnergy = sum(imEnergy_orig(:,:,1:2),3);
    elseif xys(3)>=1.7 % 1.73 is the largest we currently have in s19 and s24
       imEnergy = sum(imEnergy_orig(:,:,1:2),3);
    end

    % Average across orientations for SOC model
    imEnergyMean = zeros(size(imEnergy,1),res*res);
    for kk = 1:size(imEnergy,1)
        thisImage = imEnergy(kk,:);
        thisImage = reshape(thisImage,nrOrientations,res*res);
        imFilt1 = sum(thisImage,1); % sum across orientations
        imEnergyMean(kk,:) = imFilt1;
    end
    
    % Choose an analysis type:
%     analysisType = 'spectraRERP500';
    % analysisType = 'spectra';
%     analysisType = 'spectra500';
    analysisType = 'spectra200';

%     % load resamp_parms_even and resamp_parms_odd for even/odd:
%     dataFitName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '_evenodd.mat'];
%     load(dataFitName)
    % load resamp_parms for 100 bootstraps:
    dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)
    
    % Broadband estimate, one value per image:
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
%     ecog_bb = 10.^(resamp_parms(:,:,2)-bb_base)-1;
    ecog_bb = resamp_parms(:,:,2)-bb_base;
    % ecog_g = 10.^resamp_parms(:,:,3)-1;
    
    %% Fit SOC model on bb to confirm prf location?

    % The SOC model parameters that we will be fitting are [R C S G N] where:
    %   R is the row index of the center of the 2D Gaussian
    %   C is the column index of the center of the 2D Gaussian
    %   S is the size of the 2D Gaussian
    %   G is a gain parameter
    %   N is an exponent

    % Standard SOC parameters for visual areas:
    if v_area==1
        seed_params = [xys_pix(1) xys_pix(2) xys_pix(3)*sqrt(.18) 1 .18 .93]; 
    elseif v_area==2
        seed_params = [xys_pix(1) xys_pix(2) xys_pix(3)*sqrt(.13) 1 .13 .99]; 
    elseif v_area==3
        seed_params = [xys_pix(1) xys_pix(2) xys_pix(3)*sqrt(.12) 1 .12 .99]; 
    elseif v_area==4
        seed_params = [xys_pix(1) xys_pix(2) xys_pix(3)*sqrt(.115) 1 .115 .95]; 
    end
    
    % Get model prediction for standard parameters:
    [resultsSpace1,modelfitSpace1] = ...
        helpfit_SOC(imEnergyMean,seed_params,ecog_bb(:,1),[]);
    
    % Estimate gain for all bootstraps:
    for kk = 1:size(ecog_bb,2) % nr of boots
        % estimate gain
        bb_gain = regress(ecog_bb(:,kk),modelfitSpace1);
        
        % put parameters together
        boot_SOCparams(kk,:) = [seed_params(1) seed_params(2) seed_params(3) bb_gain seed_params(5) seed_params(6)];
        
        % calculate modelfit
        [~,modelfitSpace] = ...
            helpfit_SOC(imEnergyMean,boot_SOCparams(kk,:),[],[]);
        boot_SOCfit(kk,:) = modelfitSpace;
        
        % train performance
        boot_SOCR(kk,:) = calccod(modelfitSpace,ecog_bb(:,kk));
    end

    save(fullfile(dataDir,'soc_bids','derivatives','model_output','deriveSOCbb',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_deriveSOCbblogpower']),...
        'resultsSpace1','modelfitSpace1','xys_pix','seed_params',...
        'boot_SOCparams','boot_SOCR','boot_SOCfit')
    
end

%% 
%% Display results for one electrode

%%%%% Pick a subject:
subjects = [19,23,24];
s = 1; subj = subjects(s);

% %%%% Pick an electrode:
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3

analysisType = 'spectra200';
% modelType = 'deriveSOCbblogpower';
modelType = 'deriveSOCbbpower';

elec = 109;
res = 240;

% load model fit
load(fullfile(dataDir,'soc_bids','derivatives','model_output','deriveSOCbb',...
    ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]))

% load ecog data:
dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
    'ses-01','derivatives','ieeg',...
    ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
load(dataFitName)
    
% ecog power
if isequal(modelType,'deriveSOCbbpower')
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = mean(10.^(resamp_parms(:,:,2)-bb_base)-1,2);
    ecog_g = 10.^resamp_parms(:,:,3)-1;
elseif isequal(modelType,'deriveSOCbblogpower')
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = mean(resamp_parms(:,:,2)-bb_base,2);
    ecog_g = mean(resamp_parms(:,:,3),2);
end
figure('Position',[0 0 1200 160])

ylims = [[min(ecog_bb(:)) max(ecog_bb(:))];...
         [0 max([ecog_g(:); boot_SOCfit(:)])]];

subplot(1,3,1:2),hold on
bar(ecog_bb,1,'b','EdgeColor',[0 0 0]);
plot(boot_SOCfit' ,'g')
ylim(ylims(1,:))
title(['elec ' int2str(elec)])
% plot stimulus cutoffs
stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
for k=1:length(stim_change)
    plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
end
xlim([0 87])
set(gca,'YTick',[0:.2:1])
ylabel('bb')

subplot(1,3,3)
%%% LOOK AT WHERE THE GAUSSIAN IS
[~,xx,yy] = makegaussian2d(res,2,2,2,2);
% gau = makegaussian2d(res,results.params(1),results.params(2),results.params(3)/sqrt(results.params(5)),results.params(3)/sqrt(results.params(5)),xx,yy,0,0);
gau = makegaussian2d(res,xys_pix(1),xys_pix(2),xys_pix(3),xys_pix(3),xx,yy,0,0);
imagesc(gau);
colorbar
axis image
hold on
% look at the prf from the SOC fit:
numPoints = 50;
c.th = linspace(0,2*pi, numPoints);
for kk = 1:size(boot_SOCparams,1)
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*boot_SOCparams(kk,3)./sqrt(boot_SOCparams(kk,5)));
    plot(c.x + boot_SOCparams(kk,2), c.y + boot_SOCparams(kk,1), 'k') % this is just reversed because plot and imagesc are opposite, checked this with contour
end
title('bar pRF (color) and SOC pRF (black)')

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','model_output','deriveSOCbb',...
        ['sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_' modelType]))
print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','model_output','deriveSOCbb',...
        ['sub-' int2str(subj) '_' analysisType '_el' int2str(elec) '_' modelType]))


