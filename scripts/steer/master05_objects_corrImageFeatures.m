
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
load(['./data/task-objects_stimuli-repeat_steerFilt01.mat'],'imEnergy')
% load(['./data/task-objects_stimuli-nonrepeat_steerFilt01.mat'],'imEnergy')
% absolute contrast deviding my mean across all images:
% load(['./data/task-objects_stimuli-repeat_steerFilt02.mat'],'imEnergy')

% % now to view an image, one can use:
% im_nr = 1;
% a = reshape(imEnergy(im_nr,:),sqrt(size(imEnergy,2)/8),sqrt(size(imEnergy,2)/8),8);
% figure,imagesc(a(:,:,4))


%%
%% Reliability different ECoG analyses
%%

electrodes = [66 67 68 69];
    
analysisTypes = {'spectra200','spectra','spectra500'}; % window length for epoch 0-1000 ms: 200, 300 and 500

cod_bb = zeros(length(electrodes),length(analysisTypes));
cod_g = zeros(length(electrodes),length(analysisTypes));
cod_a = zeros(length(electrodes),length(analysisTypes));

for elInd = 1:length(electrodes)
    elec = electrodes(elInd);
    
    for kk = 1:length(analysisTypes)

        %%%% Choose an analysis type for Broadband:
        analysisType = analysisTypes{kk};
        %%%% Load the fitting results:
        dataFitName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
        load(dataFitName)
        resamp_parms = resamp_parms_allrepeat;
        %%%% Broadband estimate, one value per image:
        bb_base = resamp_parms(1,6); % from the baseline is the same for resamp_parms(:,6)
    %     ecog_bb = squeeze(resamp_parms(:,2))-bb_base;
    %     ecog_g = squeeze(resamp_parms(:,3));
    %     ecog_a = squeeze(resamp_parms(:,7));

        cod_bb(elInd,kk) = calccod(resamp_parms_evenrepeat(:,2)-bb_base,resamp_parms_oddrepeat(:,2)-bb_base);
        cod_g(elInd,kk) = calccod(resamp_parms_evenrepeat(:,3),resamp_parms_oddrepeat(:,3));
        cod_a(elInd,kk) = calccod(resamp_parms_evenrepeat(:,7),resamp_parms_oddrepeat(:,7));
    end
end

figure('Position',[0 0 250 400])

subplot(2,1,1),hold on
title('Reliability for different window length')
bar(cod_bb./100)
xlabel('electrode')
ylabel('BB COD (R)')
ylim([-0.1 1]), set(gca,'XTick',[1 2 3 4],'XTickLabel',electrodes), grid on

subplot(2,1,2),hold on
bar(cod_g./100)
xlabel('electrode')
ylabel('Gamma COD (R)')
ylim([-0.1 1]), set(gca,'XTick',[1 2 3 4],'XTickLabel',electrodes), grid on
legend('200','300','500')

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/spectralReliability/sub-' subj '_CODdifferentAnalyses'])
% print('-dpng','-r300',['./figures/spectralReliability/sub-' subj '_CODdifferentAnalyses'])

figure('Position',[0 0 250 400])

subplot(2,1,1),hold on
title('Reliability for different window length')
bar(cod_a./100)
xlabel('electrode')
ylabel('Alpha COD (R)')
ylim([-0.1 1]), set(gca,'XTick',[1 2 3 4],'XTickLabel',electrodes), grid on


%%
%% get PRF info and load ECoG data from one electrode
%%

elec = 68;
   
%%%% Choose an analysis type for Broadband:
analysisType = 'spectra'; % 0-1000 ms, 300ms window: use for BB
%%%% Load the fitting results:
dataFitName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
load(dataFitName)
resamp_parms = resamp_parms_allrepeat;
% resamp_parms = resamp_parms_nonrepeat;
%%%% Broadband estimate, one value per image:
bb_base = resamp_parms(1,6); % from the baseline is the same for resamp_parms(:,6)
ecog_bb = squeeze(resamp_parms(:,2))-bb_base;
% ecog_a = squeeze(resamp_parms(:,7));
cod_bb = calccod(resamp_parms_evenrepeat(:,2)-bb_base,resamp_parms_oddrepeat(:,2)-bb_base);

%%%% Choose an analysis type for Broadband:
analysisType = 'spectra500'; % 0-1000 ms, 500ms window: use for gamma
%%%% Load the fitting results:
dataFitName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
load(dataFitName)
resamp_parms = resamp_parms_allrepeat;
% resamp_parms = resamp_parms_nonrepeat;
%%%% Gamma estimate, one value per image:
ecog_g = squeeze(resamp_parms(:,3));
% ecog_a = squeeze(resamp_parms(:,7));
cod_g = calccod(resamp_parms_evenrepeat(:,3),resamp_parms_oddrepeat(:,3));

disp(['even/odd cod for bb: ' num2str(cod_bb) ' , gamma: ' num2str(cod_g)])
% clear resamp_parm*

% get prf info
res = sqrt(size(imEnergy,2)/8);  % resolution of the pre-processed stimuli
[gau,xx,yy] = makegaussian2d(res,2,2,2,2);

%%% add pRF fit from bar-task
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

% Get the Prf center
[~,temp_y_ind] = max(sum(G_el,1));
[~,temp_x_ind] = max(sum(G_el,2));
prf_x = xx_1(temp_x_ind,temp_y_ind);
prf_y = yy_1(temp_x_ind,temp_y_ind);
prf_s = prf_info(3)*(floor((res+1)/2)/5.5); % sigma in pixels
clear temp_x_max temp_x_ind temp_y_ind % housekeeping
pp = [prf_y prf_x max([2*prf_s 10])]; 
% Make the super small prf's a little bit larger...

% Define prf center Gaussian:
gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
    pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

% Define prf surround Gaussian:
% surndIncrease = 1.5;
% gaufun2 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),surndIncrease*pp(3),...
%     surndIncrease*pp(3),xx,yy,0,0)/(2*pi*surndIncrease*pp(3)^2));
% surroundfun = @(wts2,wts1) restrictrange(bsxfun(@minus,wts2,wts1),0,Inf);

% % Define a model to calculate the contrast within the prf:
% modelfun = @(pp,dd) dd*gaufun1(pp);
% % modelfun = @(pp,dd) dd*surroundfun(gaufun2(pp),gaufun1(pp));

% Plot prf
figure,
imagesc(reshape(gaufun1(pp),res,res),[0 max(gaufun1(pp))]),hold on
contour(xx_1,yy_1,G_el,contourVal,'k')
title('prf')

% Get mean image energies across orientations:
imEnergyMean = zeros(size(imEnergy,1),res*res);
for kk = 1:size(imEnergy,1)
    thisImage = imEnergy(kk,:);
    thisImage = reshape(thisImage,res*res,8);
    imFilt1 = sum(thisImage,2); % Sum across orientations
    imEnergyMean(kk,:) = imFilt1;
end

%% Get image values for all images

out_raw = getRawImageVals(ims_repeat,pp,res,xx,yy);
out_ce = getCEImageVals(imEnergy,pp,res,xx,yy);

out_fourier = getFourierImageVals(ims_repeat,pp,res,xx,yy);

out_norm = getDivNormImageVals(imEnergy,pp,res,xx,yy);

%% just some tests:
allVals = [out_ce.prf1CE out_ce.prf2CE out_ce.prf3CE out_ce.prf4CE ...
    out_ce.totalCE out_ce.prfVarCE out_ce.prfSurCE out_raw.minMaxIM out_raw.varIM ...
    out_raw.CEPrfIm out_raw.varPrfIM out_raw.test out_fourier.rptIm];

% allVals = [out_ce.prf2CE...
%     out_ce.prfVarCE...
%     out_raw.varIM...
%     out_raw.test];

corr(allVals,ecog_g).^2

stats = regstats(ecog_g,allVals);
[stats.rsquare stats.adjrsquare]

allVals2 = [out_ce.prf2CE...
    out_ce.prfVarCE...
    out_raw.varIM];
stats2 = regstats(ecog_g,allVals2);
[stats2.rsquare stats2.adjrsquare]

allVals2 = [out_norm.NormR(:,2)...
    out_raw.varIM];
stats2 = regstats(ecog_g,allVals2);
[stats2.rsquare stats2.adjrsquare]


%%
%% Loop across electrodes and plot
%%

clear out_ecog out_raw out_ce out_fourier out_norm
electrodes = [66 67 68 69];

cod_bb = zeros(length(electrodes),1);
cod_g = zeros(length(electrodes),1);

for elInd = 1:length(electrodes)

    elec = electrodes(elInd);

    %%%% Choose an analysis type for Broadband:
%     analysisType = 'spectra'; % 0-1000 ms, 300ms window: use for BB
    analysisType = 'spectra100'; % 0-1000 ms, 300ms window: use for BB
    %%%% Load the fitting results:
    dataFitName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
    load(dataFitName)
    %%%% Broadband estimate, one value per image:
    bb_base = resamp_parms_oddrepeat(1,6); % from the baseline is the same for resamp_parms(:,6)
    ecog_bb = squeeze(resamp_parms_oddrepeat(:,2))-bb_base;
    ecog_bb_test = squeeze(resamp_parms_evenrepeat(:,2))-bb_base;
    % ecog_a = squeeze(resamp_parms(:,7));
    cod_bb(elInd) = calccod(ecog_bb,ecog_bb_test);

    %%%% Choose an analysis type for Broadband:
%     analysisType = 'spectra500'; % 0-1000 ms, 500ms window: use for gamma
    analysisType = 'spectra100'; % 0-1000 ms, 500ms window: use for gamma
    %%%% Load the fitting results:
    dataFitName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
    load(dataFitName)
    %%%% Gamma estimate, one value per image:
    ecog_g = squeeze(resamp_parms_oddrepeat(:,3));
    ecog_g_test = squeeze(resamp_parms_evenrepeat(:,3));
    % ecog_a = squeeze(resamp_parms(:,7));
    cod_g(elInd) = calccod(ecog_g,ecog_g_test); % x,y

    disp(['even/odd cod for bb: ' num2str(cod_bb(elInd)) ' , gamma: ' num2str(cod_g(elInd))])
    % clear resamp_parm*
    out_ecog(elInd).ecog_g = ecog_g;
    out_ecog(elInd).ecog_g_test = ecog_g_test;
    out_ecog(elInd).ecog_bb = ecog_bb;
    out_ecog(elInd).ecog_bb_test = ecog_bb_test;

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
    out_raw(elInd) = getRawImageVals(ims_repeat,pp,res,xx,yy);
    out_ce(elInd) = getCEImageVals(imEnergy,pp,res,xx,yy);
    out_fourier(elInd) = getFourierImageVals(ims_repeat,pp,res,xx,yy);
    out_norm(elInd) = getDivNormImageVals(imEnergy,pp,res,xx,yy);
    
end

%% Correlate gamma/bb with CE in prf

% Regression only prf CE
stats_g = zeros(length(electrodes),2); % 12 output values from raw and ce
stats_bb = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE];
       
    allVals = [...
            out_ce(elInd).prf1CE...
            out_ce(elInd).prfVarCE...
            out_raw(elInd).varIM...
            ...
            ];

%     allVals = [...
%             out_ce(elInd).prf1CE ...
%             out_ce(elInd).prfVarCE...
%             out_fourier(elInd).rptPrf...
%             out_norm(elInd).NBF...
%             ];

% %     Regress and r2
%     thisStats = regstats(out_ecog(elInd).ecog_g,allVals);
%     stats_g(elInd,:) = [thisStats.rsquare thisStats.adjrsquare];
%     thisStats = regstats(out_ecog(elInd).ecog_bb,allVals);
%     stats_bb(elInd,:) = [thisStats.rsquare thisStats.adjrsquare];
    
    % Regress and COD with odd
    thisStats = regstats(out_ecog(elInd).ecog_g,allVals);
    stats_g(elInd,:) = calccod(thisStats.yhat,out_ecog(elInd).ecog_g_test)./100;

    thisStats = regstats(out_ecog(elInd).ecog_bb,allVals);
    stats_bb(elInd,:) = calccod(thisStats.yhat,out_ecog(elInd).ecog_bb_test)./100;
end

%%% TODO: change figure to COD
figure('Position',[0 0 300 250])
subplot(1,2,1),hold on
bar(stats_g(:,1))
plot(restrictrange(cod_g./100,0,1),'r*')
xlim([0 5]),ylim([0 1]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.1:1]), grid on
xlabel('electrode')
ylabel('COD(R) - CE in prf')
title('gamma')

subplot(1,2,2),hold on
bar(stats_bb(:,1))
plot(restrictrange(cod_bb./100,0,1),'r*')
xlim([0 5]),ylim([0 1]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.1:1]), grid on
title('broadband')
xlabel('electrode')

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/corr_CE_ECoG/sub-' subj '_corrCEPrfVarCEVarIm_bbgg'])
% print('-dpng','-r300',['./figures/corr_CE_ECoG/sub-' subj '_corrCEPrfVarCEVarIm_bbgg'])

%% Correlate gamma and make figures

out = zeros(length(electrodes),19); % N output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prf2CE out_ce(elInd).prf3CE out_ce(elInd).prf4CE ...%1:4
            out_ce(elInd).totalCE out_ce(elInd).prfVarCE out_ce(elInd).prfSurCE ...%5:7
            out_raw(elInd).minMaxIM out_raw(elInd).varIM ...%8:9
            out_raw(elInd).CEPrfIm out_raw(elInd).varPrfIM out_raw(elInd).test...%10:12
            out_fourier(elInd).rptIm out_fourier(elInd).rptPrf...%13:14
            out_norm(elInd).NormR...%15:18
            out_norm(elInd).NBF];%19
    for kk = 1:size(allVals,2)
        thisStats = regstats(out_ecog(elInd).ecog_g,allVals(:,kk));
        out(elInd,kk) = calccod(thisStats.yhat,out_ecog(elInd).ecog_g_test)./100;
    end
end

figure('Position',[0 0 800 800])
subplot(5,1,1),hold on
bar(out')
ylim([0 1]),set(gca,'YTick',[0:0.1:1]), grid on
xlim([0 size(out,2)+1]),set(gca,'XTick',[1:size(out,2)])

% Regression all values
subplot(5,2,3),hold on
stats = zeros(length(electrodes),2); 
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prf2CE out_ce(elInd).prf3CE out_ce(elInd).prf4CE ...
            out_ce(elInd).totalCE out_ce(elInd).prfVarCE out_ce(elInd).prfSurCE ...
            out_raw(elInd).minMaxIM out_raw(elInd).varIM ...
            out_raw(elInd).CEPrfIm out_raw(elInd).varPrfIM out_raw(elInd).test...
            out_fourier(elInd).rptIm out_fourier(elInd).rptPrf...
            out_norm(elInd).NormR...
            out_norm(elInd).NBF];
%     allVals = [out_ce(elInd).prfVarCE...
%             out_norm(elInd).NormR(:,2)];
%     allVals = [out_ce(elInd).prf1CE ...
%             out_ce(elInd).prfSurCE...
%             out_ce(elInd).prfVarCE...
%             out_raw(elInd).varIM ...
%             ...
%             out_fourier(elInd).rptPrf...
%             ];
    thisStats = regstats(out_ecog(elInd).ecog_g,allVals);
    stats(elInd,1) = calccod(thisStats.yhat,out_ecog(elInd).ecog_g_test)./100;
end
bar(stats(:,1))
ylim([0 .6]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.1:1]), grid on
title('COD(R) - all inputs')

% Regression only prf CE
subplot(5,2,4),hold on
stats = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE];
    thisStats = regstats(out_ecog(elInd).ecog_g,allVals);
    stats(elInd,1) = calccod(thisStats.yhat,out_ecog(elInd).ecog_g_test)./100;
end
bar(stats(:,1))
ylim([0 .6]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.1:1]), grid on
title('R^2 - CE in prf')

% Regression prf CE + prf variance
subplot(5,2,5),hold on
stats = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prfVarCE out_raw(elInd).varIM];
    thisStats = regstats(out_ecog(elInd).ecog_g,allVals);
    stats(elInd,1) = calccod(thisStats.yhat,out_ecog(elInd).ecog_g_test)./100;
end
bar(stats(:,1))
ylim([0 .6]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.1:1]), grid on
title('R^2 - CE in prf + variance in CE in prf + Var in image')

% Regression selection
subplot(5,2,6),hold on
stats = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prfVarCE out_raw(elInd).varIM out_raw(elInd).test];
    thisStats = regstats(out_ecog(elInd).ecog_g,allVals);
    stats(elInd,1) = calccod(thisStats.yhat,out_ecog(elInd).ecog_g_test)./100;
end
bar(stats(:,1))
ylim([0 .6]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.1:1]), grid on
title('R^2 - CE in prf + Var in CE + varIm  + rptIm ')

%% Correlate broadband and make figures

out = zeros(length(electrodes),14); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prf2CE out_ce(elInd).prf3CE out_ce(elInd).prf4CE ...
            out_ce(elInd).totalCE out_ce(elInd).prfVarCE out_ce(elInd).prfSurCE ...
            out_raw(elInd).minMaxIM out_raw(elInd).varIM ...
            out_raw(elInd).CEPrfIm out_raw(elInd).varPrfIM out_raw(elInd).test ...
            out_fourier(elInd).rptIm out_fourier(elInd).rptPrf];
    out(elInd,:) = corr(allVals,out_ecog(elInd).ecog_bb).^2;
end

figure
subplot(5,1,1),hold on
bar(out')
ylim([0 1]),set(gca,'YTick',[0:0.2:1]), grid on
xlim([0 15]),set(gca,'XTick',[1:14])

% Regression all values
subplot(5,1,2),hold on
stats = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prf2CE out_ce(elInd).prf3CE out_ce(elInd).prf4CE ...
            out_ce(elInd).totalCE out_ce(elInd).prfVarCE out_ce(elInd).prfSurCE ...
            out_raw(elInd).minMaxIM out_raw(elInd).varIM ...
            out_raw(elInd).CEPrfIm out_raw(elInd).varPrfIM out_raw(elInd).test];
    thisStats = regstats(out_ecog(elInd).ecog_bb,allVals);
    stats(elInd,:) = [thisStats.rsquare thisStats.adjrsquare];
end
bar(stats)
ylim([0 1]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.2:1]), grid on

% Regression only prf CE
subplot(5,1,3),hold on
stats = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE];
    thisStats = regstats(out_ecog(elInd).ecog_bb,allVals);
    stats(elInd,:) = [thisStats.rsquare thisStats.adjrsquare];
end
bar(stats)
ylim([0 1]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.2:1]), grid on

% Regression prf CE + prf variance
subplot(5,1,4),hold on
stats = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prfVarCE];
    thisStats = regstats(out_ecog(elInd).ecog_bb,allVals);
    stats(elInd,:) = [thisStats.rsquare thisStats.adjrsquare];
end
bar(stats)
ylim([0 1]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.2:1]), grid on

% Regression selection
subplot(5,1,5),hold on
stats = zeros(length(electrodes),2); % 12 output values from raw and ce
for elInd = 1:length(electrodes)
    allVals = [out_ce(elInd).prf1CE out_ce(elInd).prfVarCE out_raw(elInd).varIM out_raw(elInd).test];
    thisStats = regstats(out_ecog(elInd).ecog_bb,allVals);
    stats(elInd,:) = [thisStats.rsquare thisStats.adjrsquare];
end
bar(stats)
ylim([0 1]), set(gca,'XTick',[1 2 3 4]),set(gca,'YTick',[0:0.2:1]), grid on


%%
%% Test which variable adds explained variance on top of contrast energy in prf
%%

r2_start = zeros(length(electrodes),1);
r2_other = zeros(length(electrodes),18);

for elInd = 1:length(electrodes)
    start_val = [out_ce(elInd).prf1CE];
    baseStats = regstats(out_ecog(elInd).ecog_g,start_val);
    r2_start(elInd) = baseStats.adjrsquare;
    
    other_val = [out_ce(elInd).prf2CE out_ce(elInd).prf3CE out_ce(elInd).prf4CE ...% 1:3
            out_ce(elInd).totalCE out_ce(elInd).prfVarCE out_ce(elInd).prfSurCE ...% 4:6
            out_raw(elInd).minMaxIM out_raw(elInd).varIM ...% 7:8
            out_raw(elInd).CEPrfIm out_raw(elInd).varPrfIM out_raw(elInd).test...% 9 10 11
            out_fourier(elInd).rptIm out_fourier(elInd).rptPrf...% 12 13
            out_norm(elInd).NormR...% 14:17
            out_norm(elInd).NBF];% 18
        
    for kk = 1:size(other_val,2)
        thisStats = regstats(out_ecog(elInd).ecog_g,[start_val other_val(:,kk)]);
        r2_other(elInd,kk) = thisStats.adjrsquare;
    end
end

figure
subplot(2,2,1)
bar(r2_start)
title('Predictor 1 prf CE')

subplot(2,2,2),hold on
plot([1:18],zeros(1,18),'k')
% bar(r2_other-r2_start)
% plot(r2_other-r2_start)
plot((r2_other(1:3,:)-r2_start(1:3,:))')
title('added adjusted R2')
 
% r2_start = zeros(length(electrodes),1);
% r2_other = zeros(length(electrodes),9);
% 
% for elInd = 1:length(electrodes)
%     start_val = [out_ce(elInd).prf1CE out_ce(elInd).prfVarCE];
%     baseStats = regstats(out_ecog(elInd).ecog_g,start_val);
%     r2_start(elInd) = baseStats.adjrsquare;
%    
%     other_val = [out_ce(elInd).prf2CE out_ce(elInd).prf3CE out_ce(elInd).prf4CE ...
%             out_ce(elInd).totalCE out_ce(elInd).prfSurCE ...
%             out_raw(elInd).minMaxIM out_raw(elInd).varIM ...
%             out_raw(elInd).CEPrfIm out_raw(elInd).varPrfIM out_raw(elInd).test];
%         
%     for kk = 1:size(other_val,2)
%         thisStats = regstats(out_ecog(elInd).ecog_g,[start_val other_val(:,kk)]);
%         r2_other(elInd,kk) = thisStats.adjrsquare;
%     end
% end
% subplot(2,2,3)
% bar(r2_other-r2_start)
% title('added adjusted R2')

