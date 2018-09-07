%
% Correlate image contrast of patterns filtered with gabor filters with
% gamma/broadband power.
%
% Dora Hermes, 2017

clear all
% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath('/Users/dora/Documents/git/ecogBasicCode/')
addpath(genpath('/Users/dora/Documents/m-files/knkutils'));

%% Load preprocessed images and divisive normalization:

load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

% compute the square root of the sum of the squares of the outputs of
% quadrature-phase filter pairs (this is the standard complex-cell energy model).
% after this step, stimulus is images x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2));

% % compute the population term in the divisive-normalization equation.
% % this term is simply the average across the complex-cell outputs
% % at each position (averaging across orientation).
% stimulusPOP = blob(stimulus,2,8)/8;
% 
% % repeat the population term for each of the orientations
% stimulusPOP = upsamplematrix(stimulusPOP,8,2,[],'nearest');
% 
% % apply divisive normalization to the complex-cell outputs.  there are two parameters
% % that influence this operation: an exponent term (r) and a semi-saturation term (s).
% % the parameter values specified here were determined through a separate fitting
% % procedure (see paper for details).  for the purposes of this script, we will
% % simply hard-code the parameter values here and not worry about attempting to fit
% % the parameters.
% r = 1;
% s = 0.5;
% stimulus = stimulus.^r ./ (s.^r + stimulusPOP.^r);
% clear stimulusPOP;
% 
% % sum across orientation.  after this step, stimulus is images x positions.
% imEnergyMean = blob(stimulus,2,8);


%% Load ECoG data and fit

%%%%% Pick a subject:
subjects = [19,23,24];
s = 3; subj = subjects(s);

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
% electrodes = [109];

res = sqrt(size(stimulus,2)/nrOrientations);  % resolution of the pre-processed stimuli

for el = 1:length(electrodes)
    elec = electrodes(el);

    [v_area,xys,roi_labels] = subj_prf_info(subj,elec);
    % Convert xys from degrees to pixels
    xys_pix = [res./2 res./2 0] + res.*(xys./im_deg);
    xys_pix(1:2) = [res-xys_pix(2) xys_pix(1)]; % make sure that it matches images 
    
    [prf_s] = xy2prfsize(xys,v_area);
    prf_s_pix = res.*(prf_s./im_deg);
    
    % Choose an analysis type:
%     analysisType = 'spectraRERP500';
    % analysisType = 'spectra';
%     analysisType = 'spectra500';
    analysisType = 'spectra200';
    
    % load ECoG data:
    dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)
    
    % Broadband power estimate, one value per image:
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = 10.^(resamp_parms(:,:,2)-bb_base)-1;
    ecog_g = 10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1;
    % ecog_g = resamp_parms(:,:,3)./resamp_parms(:,:,5);
    
    %% Load SOC model results to get x and y and sigma/sqrt(n)
    
    modelType = 'fitSOCbbpower2';
    load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','fitSOCbb',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'xys_pix','seeds',...
        'cross_SOCparams','cross_SOCestimate')
    % take the median of the fitting parameters, is this good?
    cross_SOCparams(cross_SOCparams(:,6)<0,6) = 0; % restrictrange c min at 0
    cross_SOCparams(cross_SOCparams(:,6)>1,6) = 1; % restrictrange c max at 1
    cross_SOCparams(:,3) = abs(cross_SOCparams(:,3)); % size>0
    medianParams = median(cross_SOCparams,1);
    
    %% Fit gamma model.

    % run through sizes and fit gains:
    prf_sizes = 1;%.5:.5:5; % factor to multiply prf size

    % Initalize outputs:
    cross_NBFparams = zeros(size(ecog_bb,1),4,length(prf_sizes));
    cross_NBFestimate = zeros(size(ecog_bb,1),1,length(prf_sizes));
    
    % Estimate gain, leave one out every time for cross validation
    % seed_params = [medianParams(1:2) medianParams(3)./sqrt(medianParams(5)) 1];

    seed_params = [medianParams(1:2) prf_s_pix 1];    

    for kk = 1:size(ecog_bb,1) % number of stimuli, leave out kk   
        % training stimuli (kk is left out)
        trainSet = setdiff([1:size(ecog_bb,1)],kk);
        
        % run model through training stimuli
        [~,modelfit] = helpfit_NBFsimple(stimulus(trainSet,:),seed_params,[],[]);
        % get the gain
        B = regress(mean(ecog_g(trainSet,:),2),modelfit);
        
        % put estimated model parameters in matrix to save
        cross_NBFparams(kk,:,1) = [seed_params(1:3) B];
        
        % run through left out data with specific gain to get fit for
        % leftout stimulus
        [~,kkEstimate] = helpfit_NBFsimple(stimulus(kk,:),[seed_params(1:3) B],[],[]);
        cross_NBFestimate(kk,1,1) = kkEstimate;
        
    end

    save(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveNBF',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_NBFsimple2']),...
        'xys_pix','seed_params','prf_sizes',...
        'cross_NBFparams','cross_NBFestimate')
    % NBFsimple uses sigma./sqrt(n) as seed for size
    % NBFsimple2 uses derived prf size from eccentricity as seed for size
    % NBFsimple3 uses derived prf size and log power
    
end


%% Display model results for one electrode

%%%%% Pick a subject:
subjects = [19,23,24];
s = 3; subj = subjects(s);

% %%%% Pick an electrode:
% electrodes = [107 108 109 115 120 121]; % S1
% electrodes = [53 54]; % S2
% electrodes = [45 46]; % S3

analysisType = 'spectra200';
modelType = 'NBFsimple2';

elec = 46;
res = 240;

% load model fit
load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveNBF',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'xys_pix','seed_params','prf_sizes',...
        'cross_NBFparams','cross_NBFestimate')

% load ecog data:
dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
    'ses-01','derivatives','ieeg',...
    ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
load(dataFitName)

% ecog power
bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
ecog_bb = mean(10.^(resamp_parms(:,:,2)-bb_base)-1,2);
ecog_g = mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2);

ylims = [min(ecog_g(:)) max([ecog_g(:); cross_NBFestimate(:)+.1])];

figure('Position',[0 0 1200 160])

subplot(1,3,1:2),hold on
bar(ecog_g,1,'b','EdgeColor',[0 0 0]);
plot(cross_NBFestimate' ,'g','LineWidth',2)
% plot stimulus cutoffs
stim_change = [38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
for k = 1:length(stim_change)
    plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
end
xlim([0 87]), ylim(ylims(1,:))
set(gca,'XTick',[1:86],'XTickLabel',[],'YTick',[0:ceil(max(ecog_g(:))/4):max(ecog_g(:))])
ylabel('gamma')

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/fit' modelType  '/sub-' int2str(subj) '_' analysisType '_el' int2str(elec)])
% print('-dpng','-r300',['./figures/fit' modelType  '/sub-' int2str(subj) '_' analysisType '_el' int2str(elec)])


%% Display model results for all electrodes

%%%%% Pick a subject:
subject_ind = [19 19  19  19  19  19  24 24];
electrodes = [107 108 109 115 120 121 45 46];

figure('Position',[0 0 1400 800])
    
for ll = 1:length(electrodes)
    
    subj = subject_ind(ll);
    elec = electrodes(ll);

    analysisType = 'spectra200';
    modelType = 'NBFsimple2';

    res = 240;

    % load model fit
    load(fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveNBF',...
            ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
            'xys_pix','seed_params','prf_sizes',...
            'cross_NBFparams','cross_NBFestimate')

    % load ecog data:
    dataFitName = fullfile(dataDir,'soc_bids',['sub-' int2str(subj)],...
        'ses-01','derivatives','ieeg',...
        ['sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '.mat']);
    load(dataFitName)

    % ecog power
    bb_base = resamp_parms(1,1,6); % from the baseline is the same for resamp_parms(:,:,6)
    ecog_bb = mean(10.^(resamp_parms(:,:,2)-bb_base)-1,2);
    ecog_bb_yneg = median(10.^(resamp_parms(:,:,2)-bb_base)-1,2)-...
        quantile(10.^(resamp_parms(:,:,2)-bb_base)-1,.16,2);
    ecog_bb_ypos = quantile(10.^(resamp_parms(:,:,2)-bb_base)-1,.84,2)-...
        median(10.^(resamp_parms(:,:,2)-bb_base)-1,2);

    ecog_g = mean(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2);
    ecog_g_yneg = median(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2)-...
        quantile(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,.16,2);
    ecog_g_ypos = quantile(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,.84,2)-...
        median(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5))-1,2);
%     ecog_g = mean((resamp_parms(:,:,3)./resamp_parms(:,:,5)),2);
%     ecog_g_yneg = median((resamp_parms(:,:,3)./resamp_parms(:,:,5)),2)-...
%         quantile((resamp_parms(:,:,3)./resamp_parms(:,:,5)),.16,2);
%     ecog_g_ypos = quantile(10.^(resamp_parms(:,:,3)./resamp_parms(:,:,5)),.84,2)-...
%         median((resamp_parms(:,:,3)./resamp_parms(:,:,5)),2);


    ylims = [min(ecog_g(:)) max([ecog_g(:); cross_NBFestimate(:)+.1])];

    subplot(8,2,2*ll-1),hold on

    bar(ecog_g,1,'b','EdgeColor',[0 0 0]);
    plot(cross_NBFestimate' ,'g','LineWidth',2)
    % plot stimulus cutoffs
    stim_change = [38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
    for k = 1:length(stim_change)
        plot([stim_change(k) stim_change(k)],ylims(1,:),'Color',[.5 .5 .5],'LineWidth',2)
    end
    xlim([0 87]), ylim(ylims(1,:))
    set(gca,'XTick',[1:86],'XTickLabel',[],'YTick',[0:ceil(max(ecog_g(:))/4):max(ecog_g(:))])
    ylabel('gamma')

end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveNBF',...
        [analysisType '_allel_' modelType]))
print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt','deriveNBF',...
        [analysisType '_allel_' modelType]))

    
    
%%
%%
%% OLD code
%%
%%
%%
figure('Position',[0 0 1000 200]),
stims_plot = 74:78;
for s = 1:length(stims_plot)
    subplot(1,length(stims_plot),s)
    imagesc(reshape(imEnergyMean(stims_plot(s),:),res,res))
    hold on
    plot(results.params(2),results.params(1),'k.','MarkerSize',10)
    title(['stim ' int2str(stims_plot(s)) ' prf location el ' int2str(elec)])
    axis image
end

% set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/sub-' int2str(subj) '_exampleSpace_el' int2str(elec)])
% print('-dpng','-r300',['./figures/sub-' int2str(subj) '_exampleSpace_el' int2str(elec)]


%% Display NBF model results for all electrodes

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
modelType = 'NBFsimple';
res = 240;
nrOrientations = 8;

figure('Position',[0 0 1200 800])
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
        
        % load model fit
        load(['./data/model_output/sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType])       
        cm = jet(size(prf_sizes,2));
        
        % load ecog data:
        dataFitName = [dataRootPath '/sub-' int2str(subj) '/ses-01/derivatives/ieeg/sub-' int2str(subj) '_task-soc_allruns_' analysisType '_fitEl' int2str(elec) '_evenodd.mat'];
        load(dataFitName)
        % Test models on odd
        resamp_parms = resamp_parms_odd;
        bb_base = resamp_parms(1,6); % from the baseline is the same for resamp_parms(:,:,6)
        ecog_bb = resamp_parms(:,2)-bb_base;
        ecog_g = resamp_parms(:,3);

        subplot(nr_els,4,4*el_nr-3:4*el_nr-2),hold on % prf locations 
        
        bar(ecog_g,1,'FaceColor',[.5 .5 .5],'EdgeColor',[0 0 0]);
        
        % all fits:
        for kk = 1:size(prf_sizes,2)
            plot(NBF_modelfit(kk,:),'Color',cm(kk,:),'LineWidth',1) 
        end
        ylim([min([ecog_g(:); NBF_modelfit(:)]) max([ecog_g(:); NBF_modelfit(:)])]);
        
        % plot stimulus cutoffs
        stim_change=[38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
        for k=1:length(stim_change)
            plot([stim_change(k) stim_change(k)],[min([ecog_g(:); NBF_modelfit(:)]) max([ecog_g(:); NBF_modelfit(:)])],...
                'Color',[0 0 0],'LineWidth',2)
        end
        xlim([0 87])
        set(gca,'YTick',[0:.2:1],'XTick',[])
        ylabel(['g el ' int2str(elec)])
    end
end

% prf locations 

el_nr = 0;
for s = 1:length(s_nrs)
    subj = subjects(s_nrs(s));

    for el = 1:length(electrodes{s_nrs(s)})
        el_nr = el_nr+1;
        elec = electrodes{s_nrs(s)}(el);

        % load model fit
        load(['./data/model_output/sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType])

        subplot(nr_els,4,4*el_nr-1),hold on % prf locations 
        plot([res/2 res/2],[0 res],'k')
        plot([0 res],[res/2 res/2],'k')

        % look at the prf from the SOC fit:
        numPoints = 50;
        c.th = linspace(0,2*pi, numPoints);
        for kk = 1:size(NBF_params,1)
            [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*NBF_params(kk,3));
            plot(c.x + NBF_params(kk,2), c.y + NBF_params(kk,1),...
                'Color',cm(kk,:)) % this is just reversed because plot and imagesc are opposite, checked this with contour
        end

        xlim([0 res+1]),ylim([0 res+1])
%         set(gca,'XTick',[.25*res .5*res .75*res],'XTickLabel',{'10','0','10'},...
%             'YTick',[.25*res .5*res .75*res],'YTickLabel',{'10','0','10'},...
%             'Ydir','reverse')
        set(gca,'XTick',[.25*res .5*res .75*res],'XTickLabel',[],...
            'YTick',[.25*res .5*res .75*res],'YTickLabel',[],...
            'Ydir','reverse')
        axis square
    end
end
 
% parameters 
r_crossval = zeros(nr_els,size(prf_sizes,2),3);

el_nr = 0;
for s = 1:length(s_nrs)
    subj = subjects(s_nrs(s));

    for el = 1:length(electrodes{s_nrs(s)})
        el_nr = el_nr+1;
        elec = electrodes{s_nrs(s)}(el);

        % load model fit
        load(['./data/model_output/sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType])
        
        r_crossval(el_nr,:,:) = NBF_R;
        
        subplot(nr_els,4,4*el_nr),hold on % prf locations 
        plot(prf_sizes,NBF_R(:,1)./100,'Color',[.7 .7 .7])
        for kk = 1:size(prf_sizes,2)
%             bar(prf_sizes(kk),NBF_R(kk,1),'FaceColor',cm(kk,:))
            plot(prf_sizes(kk),NBF_R(kk,1)./100,'.','Color',cm(kk,:))
        end
        
        plot(prf_sizes,NBF_R(:,2)./100,'Color',[.7 .7 .7])
        for kk = 1:size(prf_sizes,2)
%             bar(prf_sizes(kk),NBF_R(kk,2),'FaceColor',cm(kk,:))
            plot(prf_sizes(kk),NBF_R(kk,2)./100,'*','Color',cm(kk,:))
        end       

        plot(prf_sizes,NBF_R(:,3),'Color',[.7 .7 .7])
        for kk = 1:size(prf_sizes,2)
%             bar(prf_sizes(kk),NBF_R(kk,3),'FaceColor',cm(kk,:))
            plot(prf_sizes(kk),NBF_R(kk,3),'o','Color',cm(kk,:))
        end
        ylim([0 1])
        xlim([0 max(prf_sizes)+.5])
        set(gca,'XTick',prf_sizes)
        ylabel('R or r^2')

    end
end

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['./figures/fit' modelType  '/allsub_' analysisType])
print('-dpng','-r300',['./figures/fit' modelType  '/allsub_' analysisType])

%%
% r_crossval
figure('Position',[0 0 250 500])
subplot(3,1,1),hold on
for kk = 1:size(prf_sizes,2)
    bar(prf_sizes(kk),squeeze(mean(r_crossval(:,kk,1),1)),.4,'FaceColor',[1 1 1],'EdgeColor',cm(kk,:))
    plot(prf_sizes(kk),squeeze(r_crossval(:,kk,1)),'.','Color',cm(kk,:))
end
ylabel('R')
ylim([0 100]),xlim([0 max(prf_sizes)+.5])
set(gca,'XTick',prf_sizes)

subplot(3,1,2),hold on
for kk = 1:size(prf_sizes,2)
    bar(prf_sizes(kk),squeeze(mean(r_crossval(:,kk,2),1)),.4,'FaceColor',[1 1 1],'EdgeColor',cm(kk,:))
    plot(prf_sizes(kk),squeeze(r_crossval(:,kk,2)),'.','Color',cm(kk,:))
end
ylabel('R, mean subtracted')
ylim([0 100]),xlim([0 max(prf_sizes)+.5])
set(gca,'XTick',prf_sizes)

subplot(3,1,3),hold on
for kk = 1:size(prf_sizes,2)
    bar(prf_sizes(kk),squeeze(mean(r_crossval(:,kk,3),1)),.4,'FaceColor',[1 1 1],'EdgeColor',cm(kk,:))
    plot(prf_sizes(kk),squeeze(r_crossval(:,kk,3)),'.','Color',cm(kk,:))
end
ylabel('Pearson r^2')
ylim([0 1]),xlim([0 max(prf_sizes)+.5])
set(gca,'XTick',prf_sizes)
xlabel('pRF size (x SOC pRF fit)')

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['./figures/fit' modelType  '/allsub_' analysisType '_Rs'])
print('-dpng','-r300',['./figures/fit' modelType  '/allsub_' analysisType '_Rs'])

