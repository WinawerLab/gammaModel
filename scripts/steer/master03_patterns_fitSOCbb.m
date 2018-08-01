
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
    ecog_bb = 10.^(resamp_parms(:,:,2)-bb_base)-1;
    % ecog_g = 10.^resamp_parms(:,:,3)-1;
    
    %% Fit SOC model on bb to confirm prf location?

    % The SOC model parameters that we will be fitting are [R C S G N] where:
    %   R is the row index of the center of the 2D Gaussian
    %   C is the column index of the center of the 2D Gaussian
    %   S is the size of the 2D Gaussian
    %   G is a gain parameter
    %   N is an exponent

    seeds = [(1+res)/2 (1+res)/2 res/4*sqrt(0.5) 1 .5 .5];
    [resultsSpace1,modelfitSpace1] = helpfit_SOC(imEnergyMean([1:38 59:68 74:78],:),[],mean(ecog_bb([1:38 59:68 74:78],:),2),seeds);
    
    seeds = [resultsSpace1.params(1:4) .5 .5];
    boot_SOCparams = zeros(size(ecog_bb,2),size(seeds,2));
    boot_SOCR = zeros(size(ecog_bb,2),1);
    boot_SOCfit = zeros(size(ecog_bb,2),size(ecog_bb,1));
    for kk = 1:size(ecog_bb,2) % nr of boots
        % fit:
        [resultsSpace,modelfitSpace] = helpfit_SOC(imEnergyMean,[],ecog_bb(:,kk),seeds);
        
%         % get values from  and get gain:
%         [resultsSpace,modelfitSpace] = helpfit_SOC(imEnergyMean,params,[],[]);
%         
        boot_SOCparams(kk,:) = resultsSpace.params;
        boot_SOCR(kk,:) = resultsSpace.trainperformance;
        boot_SOCfit(kk,:) = modelfitSpace;
    end

%     save(['./data/model_output/sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_SOCbbpower'],...
%         'resultsSpace1','modelfitSpace1','xys_pix','seeds',...
%         'boot_SOCparams','boot_SOCR','boot_SOCfit')
    
end

%% quick plot and test a different model

params = mean(boot_SOCparams,1);
params(3) = 10;
params(5) = 0.18;
params(6) = 0.93;
[~,modelfit] = helpfit_SOC(imEnergyMean,params,[],[]);


%% 
%% Display SOC model results for one electrode

dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
%%%%% Pick a subject:
subjects = [19,23,24];
s = 1; subj = subjects(s);

% %%%% Pick an electrode:
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3

analysisType = 'spectra200';
modelType = 'SOC';

elec = 107;
res = 240;

% load model fit
load(['./data/model_output/sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType])

% load ecog data:
dataFitName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '_evenodd.mat'];
load(dataFitName)
    
% Test models on odd
resamp_parms = resamp_parms_odd;
bb_base = resamp_parms(1,6); % from the baseline is the same for resamp_parms(:,:,6)
ecog_bb = resamp_parms(:,2)-bb_base;
ecog_g = resamp_parms(:,3);

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

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/fit' modelType  '/sub-' int2str(subj) '_' analysisType '_el' int2str(elec)])
% print('-dpng','-r300',['./figures/fit' modelType  '/sub-' int2str(subj) '_' analysisType '_el' int2str(elec)])


%% Display SOC model results for all electrodes

dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';
%%%%% Pick a subject:
subjects = [19,23,24];
s_nrs = [1 3];

% electrodes
electrodes = {[107 108 109 115 120 121],[53 54],[45 46]}; % S1, S2, S3
nr_els = 0;
for kk = s_nrs
    nr_els = nr_els + length(electrodes{kk});
end
el_colors = jet(nr_els);

analysisType = 'spectra200';
modelType = 'SOC';
res = 240;
nrOrientations = 8;

figure('Position',[0 0 1200 600])
el_nr = 0;
for s = 1:length(s_nrs)
    subj = subjects(s_nrs(s));
    if subj == 9
        im_deg = rad2deg(atan(20.7./61));
    elseif subj == 19
        im_deg = rad2deg(atan(17.9./50));
    elseif subj==23
        im_deg = rad2deg(atan(17.9./36));
    elseif subj==24
        im_deg = rad2deg(atan(17.9./45));
    end

    for el = 1:length(electrodes{s_nrs(s)})
        el_nr = el_nr+1;
        elec = electrodes{s_nrs(s)}(el);
        
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
        
        % load model fit
        load(['./data/model_output/sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType])

        % load ecog data:
        dataFitName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '_evenodd.mat'];
        load(dataFitName)
        % Test models on odd
        resamp_parms = resamp_parms_odd;
        bb_base = resamp_parms(1,6); % from the baseline is the same for resamp_parms(:,:,6)
        ecog_bb = resamp_parms(:,2)-bb_base;
        ecog_g = resamp_parms(:,3);

        subplot(nr_els,3,3*el_nr-2),hold on % prf locations 
        
        bar(ecog_bb,1,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0]);
        
        % all fits:
        plot(boot_SOCfit' ,'Color',el_colors(el_nr,:))
        ylims = [min([ecog_bb(:)]) max([ecog_bb(:)])];
        % fit for average after bootstrap - error for electrode 46 s3, because of variability in gain/size/exponent:
%         mean_params = median(boot_SOCparams,1);
%         [~,modelfit] = helpfit_SOC(imEnergyMean,mean_params);
% %         plot(modelfit ,'Color',el_colors(el_nr,:),'LineWidth',2)
%         plot(modelfit ,'k','LineWidth',2)
%         ylims = [min([ecog_bb(:); modelfit]) max([ecog_bb(:); modelfit])];
        
        ylim(ylims(1,:))
        % plot stimulus cutoffs
        stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
        for k=1:length(stim_change)
            plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
        end
        xlim([0 87])
        set(gca,'YTick',[0:.2:1],'XTick',[])
        ylabel(['bb el ' int2str(elec)])
    end
end

% prf locations 

subplot(1,3,2),hold on 
plot([res/2 res/2],[0 res],'k')
plot([0 res],[res/2 res/2],'k')

el_nr = 0;
for s = 1:length(s_nrs)
    subj = subjects(s_nrs(s));

    for el = 1:length(electrodes{s_nrs(s)})
        el_nr = el_nr+1;
        elec = electrodes{s_nrs(s)}(el);

        % load model fit
        load(['./data/model_output/sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType])

        % look at the prf from the SOC fit:
        numPoints = 50;
        c.th = linspace(0,2*pi, numPoints);
        for kk = 1:size(boot_SOCparams,1)
            [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*boot_SOCparams(kk,3)./sqrt(boot_SOCparams(kk,5)));
            plot(c.x + boot_SOCparams(kk,2), c.y + boot_SOCparams(kk,1),...
                'Color',el_colors(el_nr,:)) % this is just reversed because plot and imagesc are opposite, checked this with contour
        end
        title('SOC pRF')
    end
end
xlim([0 res+1]),ylim([0 res+1])
set(gca,'XTick',[.25*res .5*res .75*res],'XTickLabel',{'10','0','10'},...
    'YTick',[.25*res .5*res .75*res],'YTickLabel',{'10','0','10'},...
    'Ydir','reverse')
axis square

%%%% TODO: LOOP FOR ALL SUBJECTS AND PLOT CI
% parameters 

el_nr = 0;
for s = 1:length(s_nrs)
    subj = subjects(s_nrs(s));

    for el = 1:length(electrodes{s_nrs(s)})
        el_nr = el_nr+1;
        elec = electrodes{s_nrs(s)}(el);

        % load model fit
        load(['./data/model_output/sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType])
        
        subplot(3,3,3),hold on 
        bar(el_nr,mean(boot_SOCparams(:,5)),'FaceColor',[1 1 1],'EdgeColor',el_colors(el_nr,:))
        plot(el_nr,boot_SOCparams(:,5),'.','Color',el_colors(el_nr,:))
        
        subplot(3,3,6),hold on 
        bar(el_nr,mean(boot_SOCparams(:,6)),'FaceColor',[1 1 1],'EdgeColor',el_colors(el_nr,:))
        plot(el_nr,boot_SOCparams(:,6),'.','Color',el_colors(el_nr,:))
        
        subplot(3,3,9),hold on 
        bar(el_nr,mean(boot_SOCR./100),'FaceColor',[1 1 1],'EdgeColor',el_colors(el_nr,:))
        plot(el_nr,boot_SOCR./100,'.','Color',el_colors(el_nr,:))
        ylim([0 1])
    end
end

subplot(3,3,3),hold on 
title('exponent')
xlim([0 9])
set(gca,'XTick',[1:8])
subplot(3,3,6),hold on 
title('c parameter')
%%% MAKE SURE TO PLOT WITH RESTRICTRANGE
asdfqet sdf x sieb as

xlim([0 9])
set(gca,'XTick',[1:8])
subplot(3,3,9),hold on 
xlabel('electrode')
xlim([0 9])
title('R on training data (no cross validation)')
set(gca,'XTick',[1:8])

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/fit' modelType  '/allsub_' analysisType])
% print('-dpng','-r300',['./figures/fit' modelType  '/allsub_' analysisType])

