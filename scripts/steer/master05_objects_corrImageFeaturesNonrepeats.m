
% Correlate image contrast of patterns filtered with steerable pyramid
% filters with gamma/broadband power.
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

subj = '09';

dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

% Load original images:
load(['../model/data/sub-' subj '_stimuli'],'ims_repeat','ims_nonrepeat')
% figure; imagesc(makeimagestack(ims_repeat)), colormap gray

% Load image energies:
% relative contrast in each image (deviding by mean for each image):
% load(['./data/task-objects_stimuli-repeat_steerFilt01.mat'],'imEnergy')
load(['./data/task-objects_stimuli-nonrepeat_steerFilt01.mat'],'imEnergy')
% absolute contrast deviding my mean across all images:
% load(['./data/task-objects_stimuli-repeat_steerFilt02.mat'],'imEnergy')

% % now to view an image, one can use:
% im_nr = 1;
% a = reshape(imEnergy(im_nr,:),sqrt(size(imEnergy,2)/8),sqrt(size(imEnergy,2)/8),8);
% figure,imagesc(a(:,:,4))


%%
%% Loop across electrodes and plot
%%

clear out_ecog out_raw out_ce out_fourier
electrodes = [66 67 68 69];
for elInd = 1:length(electrodes)

    elec = electrodes(elInd);

    %%%% Choose an analysis type for Broadband:
    analysisType = 'spectra'; % 0-1000 ms, 300ms window: use for BB
    %%%% Load the fitting results:
    dataFitName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
    load(dataFitName)
    resamp_parms = resamp_parms_nonrepeat;
    %%%% Broadband estimate, one value per image:
    bb_base = resamp_parms(1,6); % from the baseline is the same for resamp_parms(:,6)
    ecog_bb = squeeze(resamp_parms(:,2))-bb_base;

    %%%% Choose an analysis type for Broadband:
    analysisType = 'spectra500'; % 0-1000 ms, 500ms window: use for gamma
    %%%% Load the fitting results:
    dataFitName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
    load(dataFitName)
    resamp_parms = resamp_parms_nonrepeat;
    %%%% Gamma estimate, one value per image:
    ecog_g = squeeze(resamp_parms(:,3));

    % clear resamp_parm*
    out_ecog(elInd).ecog_g = ecog_g;
    out_ecog(elInd).ecog_bb = ecog_bb;

    % Get prf info
    res = sqrt(size(imEnergy,2)/8);  % resolution of the pre-processed stimuli
    [gau,xx,yy] = makegaussian2d(res,2,2,2,2);
    %%% Get pRF fit from bar-task
    f1 = figure;
    thisFile = [dataRootPath '/sub-' subj '/ses-01/derivatives/BootstrappedPRF/sub-' subj '_Chan' int2str(elec) '_exp_bb.mat'];
    [prf_info,G,xx_1,yy_1,contourVal]=lookup_09_PRFvals(elec,thisFile);
    close(f1)
    G_el=flipud(G{1});clear G
    % convert XX and YY to current image 
    %%%% assume images without padding
    xx_1(:)=xx_1(:)*(floor((res+1)/2)/5.5); % 5.5 is 11/2 - ?images were presented at 11 degrees of visual angle?
    yy_1(:)=yy_1(:)*(floor((res+1)/2)/5.5);
    xx_1(:)=xx_1(:)+floor((res+1)/2);
    yy_1(:)=yy_1(:)+floor((res+1)/2);
    %%%% assume images with padding
    % xx_1(:)=xx_1(:)*(floor((res+1)/2)/5.9714); 
    % yy_1(:)=yy_1(:)*(floor((res+1)/2)/5.9714);
    % xx_1(:)=xx_1(:)+floor((res+1)/2);
    % yy_1(:)=yy_1(:)+floor((res+1)/2);
    % Get Prf x, y, sigma:
    [~,temp_y_ind] = max(sum(G_el,1));
    [~,temp_x_ind] = max(sum(G_el,2));
    prf_x = xx_1(temp_x_ind,temp_y_ind);
    prf_y = yy_1(temp_x_ind,temp_y_ind);
    prf_s = prf_info(3)*(floor((res+1)/2)/5.5); % sigma in pixels
    clear temp_x_max temp_x_ind temp_y_ind % housekeeping
    pp = [prf_y prf_x max([2*prf_s 10])]; % make unrealistically small prfs a bit larger

    % Get image values for all images
    out_raw(elInd) = getRawImageVals(ims_nonrepeat,pp,res,xx,yy);
    out_ce(elInd) = getCEImageVals(imEnergy,pp,res,xx,yy);
    out_fourier(elInd) = getFourierImageVals(ims_nonrepeat,pp,res,xx,yy);
    
end

%% Correlate gamma and make figures

out = zeros(length(electrodes),14); % 13 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prf2CE out_ce(elInd).prf3CE out_ce(elInd).prf4CE ...
            out_ce(elInd).totalCE out_ce(elInd).prfVarCE out_ce(elInd).prfSurCE ...
            out_raw(elInd).minMaxIM out_raw(elInd).varIM ...
            out_raw(elInd).CEPrfIm out_raw(elInd).varPrfIM out_raw(elInd).test...
            out_fourier(elInd).rptIm out_fourier(elInd).rptPrf];
    out(elInd,:) = corr(allVals,out_ecog(elInd).ecog_g).^2;
end

figure('Position',[0 0 400 800])
subplot(5,1,1),hold on
bar(out')
ylim([0 .6]),set(gca,'YTick',[0:0.1:1]), grid on
xlim([0 15]),set(gca,'XTick',[1:14])

% Regression all values
subplot(5,1,2),hold on
stats = zeros(length(electrodes),2); 
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prf2CE out_ce(elInd).prf3CE out_ce(elInd).prf4CE ...
            out_ce(elInd).totalCE out_ce(elInd).prfVarCE out_ce(elInd).prfSurCE ...
            out_raw(elInd).minMaxIM out_raw(elInd).varIM ...
            out_raw(elInd).CEPrfIm out_raw(elInd).varPrfIM out_raw(elInd).test...
            out_fourier(elInd).rptIm out_fourier(elInd).rptPrf];
    thisStats = regstats(out_ecog(elInd).ecog_g,allVals);
    stats(elInd,:) = [thisStats.rsquare thisStats.adjrsquare];
end
bar(stats)
ylim([0 .6]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.1:1]), grid on
title('R^2 - all inputs')

% Regression only prf CE
subplot(5,1,3),hold on
stats = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE];
    thisStats = regstats(out_ecog(elInd).ecog_g,allVals);
    stats(elInd,:) = [thisStats.rsquare thisStats.adjrsquare];
end
bar(stats)
ylim([0 .6]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.1:1]), grid on
title('R^2 - CE in prf')

% Regression prf CE + prf variance
subplot(5,1,4),hold on
stats = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prfVarCE out_raw(elInd).varIM];
    thisStats = regstats(out_ecog(elInd).ecog_g,allVals);
    stats(elInd,:) = [thisStats.rsquare thisStats.adjrsquare];
end
bar(stats)
ylim([0 .6]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.1:1]), grid on
title('R^2 - CE in prf + Var in CE')

% Regression selection
subplot(5,1,5),hold on
stats = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prfVarCE out_raw(elInd).varIM out_fourier(elInd).rptPrf];
    thisStats = regstats(out_ecog(elInd).ecog_g,allVals);
    stats(elInd,:) = [thisStats.rsquare thisStats.adjrsquare];
end
bar(stats)
ylim([0 .6]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.1:1]), grid on
title('R^2 - CE in prf + Var in CE + varIm  + rptIm ')