%
% Compare pRF from bar scan to the pRF from the SOC model fit with size
% derived from the eccentricity.
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

load(fullfile(dataDir,'derivatives','gaborFilt','task-soc_stimuli_gaborFilt01.mat'),'stimulus')

% compute the square root of the sum of the squares of the outputs of
% quadrature-phase filter pairs (this is the standard complex-cell energy model).
% after this step, stimulus is images x orientations*positions.
stimulus = sqrt(blob(stimulus.^2,2,2));
% the square root may not be typical for complex cell energy model


%% Load ECoG data and fit

%%%%% Pick a subject:
subjects = [19,23,24,1001];
s = 3; subj = subjects(s);

nrOrientations = 8;

if subj == 9
    im_deg = rad2deg(atan(20.7./61));
elseif subj == 19
    im_deg = rad2deg(atan(17.9./50));
    electrodes = [107 108 109 115 120 121]; % S1
elseif subj == 23
    im_deg = rad2deg(atan(17.9./36));
    % electrodes = [53 54]; % S2
elseif subj == 24
    im_deg = rad2deg(atan(17.9./45));
    electrodes = [45 46]; % S3
% elseif subj == 1001
%     im_deg = rad2deg(atan(20./60));
%     electrodes = [49 50 51 52 57 58 59 60]; % S1001 V1, V2, V3
end
% electrodes = [51];

res = sqrt(size(stimulus,2)/nrOrientations);  % resolution of the pre-processed stimuli

figure('Position',[0 0 1000 400])

for el = 1:length(electrodes)
    elec = electrodes(el);
    
    % get prf from bar task
    [v_area,roi_labels,xys] = subj_prf_info(subj,elec);
    % Convert xys from degrees to pixels:
    xys_pix = [res./2 res./2 0] + res.*(xys./im_deg);
    xys_pix(1:2) = [res-xys_pix(2) xys_pix(1)]; % make sure that it matches images 
    prf_BAR = xys_pix;
    clear xys_pix
    
    % Load SOC model results to get x and y (and sigma/sqrt(n))   
    modelType = 'fitSOCbbpower2';
    load(fullfile(dataDir,'derivatives','gaborFilt','fitSOCbb',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seeds','cross_SOCparams','cross_SOCestimate')
    % take the median of the fitting parameters, is this good?
    cross_SOCparams(cross_SOCparams(:,6)<0,6) = 0; % restrictrange c min at 0
    cross_SOCparams(cross_SOCparams(:,6)>1,6) = 1; % restrictrange c max at 1
    cross_SOCparams(:,3) = abs(cross_SOCparams(:,3)); % size>0
    medianParams = median(cross_SOCparams,1);
    
    %%%%% Derive pRF size from eccentricity %%%%%
    % Convert xys from pixels to degrees
    xys_deg = (medianParams(1:2)-[res./2 res./2])*im_deg./res;
    xys_deg(1:2) = [xys_deg(2) xys_deg(1)]; % make sure that it matches images

    % Derive pRF size from eccentricity:
    [prf_s] = xy2prfsize(xys_deg,v_area);
    prf_s_pix = res.*(prf_s./im_deg);
    
    % OV model pRF in pixels:
    prf_OVmodel = [medianParams(1:2) prf_s_pix];
    clear prf_s prf_s_pix xys_deg
    
    subplot(2,4,el)
    %%% look at where the bar task Prf Gaussian is
    [~,xx,yy] = makegaussian2d(res,2,2,2,2);
    gau = makegaussian2d(res,prf_BAR(1),prf_BAR(2),prf_BAR(3),prf_BAR(3),xx,yy,0,0);
    imagesc(gau,[0 1]);

    % plotting settings
    axis image
    hold on

    %%% look at where the SOC task Prf Gaussian is
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*prf_OVmodel(3));
    plot(c.x + prf_OVmodel(2), c.y + prf_OVmodel(1), 'k') % this is just reversed because plot and imagesc are opposite, checked this with contour
    title(['el' int2str(el)])
end


%% Load ECoG data and fit

nrOrientations = 8;
res = sqrt(size(stimulus,2)/nrOrientations);  % resolution of the pre-processed stimuli

% start the figure:
figure
%%% get the correct axes
[~,xx,yy] = makegaussian2d(res,2,2,2,2);
gau = makegaussian2d(res,res/2,res/2,res/10,res/10,xx,yy,0,0);
% plot a background
imagesc(.7+zeros(size(gau)),[0 1]);
% plotting settings
axis image
hold on
colormap gray

%%%%% Pick a subject:
subjects = [19 19 19 19 19 19 24 24]; 
electrodes = [107 108 109 115 120 121 45 46];

el_color = jet(length(electrodes));

for el = 1:length(electrodes)
    
    subj = subjects(el);
    elec = electrodes(el);
    
    % get image to degrees conversion factor
    if subj == 9
        im_deg = rad2deg(atan(20.7./61));
    elseif subj == 19
        im_deg = rad2deg(atan(17.9./50));
    elseif subj == 23
        im_deg = rad2deg(atan(17.9./36));
    elseif subj == 24
        im_deg = rad2deg(atan(17.9./45));
    elseif subj == 1001
        im_deg = rad2deg(atan(20./60));
    end
    
    % get prf from bar task
    [v_area,roi_labels,xys] = subj_prf_info(subj,elec);
    % Convert xys from degrees to pixels:
    xys_pix = [res./2 res./2 0] + res.*(xys./im_deg);
    xys_pix(1:2) = [res-xys_pix(2) xys_pix(1)]; % make sure that it matches images 
    prf_BAR = xys_pix;
    clear xys_pix
    
    % Load SOC model results to get x and y (and sigma/sqrt(n))   
    modelType = 'fitSOCbbpower2';
    load(fullfile(dataDir,'derivatives','gaborFilt','fitSOCbb',...
        ['sub' int2str(subj) '_el' int2str(elec) '_' analysisType '_' modelType]),...
        'seeds','cross_SOCparams','cross_SOCestimate')
    % take the median of the fitting parameters, is this good?
    cross_SOCparams(cross_SOCparams(:,6)<0,6) = 0; % restrictrange c min at 0
    cross_SOCparams(cross_SOCparams(:,6)>1,6) = 1; % restrictrange c max at 1
    cross_SOCparams(:,3) = abs(cross_SOCparams(:,3)); % size>0
    medianParams = median(cross_SOCparams,1);
    
    %%%%% Derive pRF size from eccentricity %%%%%
    % Convert xys from pixels to degrees
    xys_deg = (medianParams(1:2)-[res./2 res./2])*im_deg./res;
    xys_deg(1:2) = [xys_deg(2) xys_deg(1)]; % make sure that it matches images

    % Derive pRF size from eccentricity:
    [prf_s] = xy2prfsize(xys_deg,v_area);
    prf_s_pix = res.*(prf_s./im_deg);
    
    % OV model pRF in pixels:
    prf_OVmodel = [medianParams(1:2) prf_s_pix];
    clear prf_s prf_s_pix xys_deg
    
    %%% look at where the bar task Prf Gaussian is
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*prf_BAR(3));
    plot(c.x + prf_BAR(2), c.y + prf_BAR(1),'-','Color',el_color(el,:),'LineWidth',2) % this is just reversed because plot and imagesc are opposite, checked this with contour

    %%% look at where the SOC task Prf Gaussian is
    numPoints = 50;
    c.th = linspace(0,2*pi, numPoints);
    [c.x, c.y] = pol2cart(c.th, ones(1,numPoints)*prf_OVmodel(3));
    plot(c.x + prf_OVmodel(2), c.y + prf_OVmodel(1),':','Color',el_color(el,:),'LineWidth',2) % this is just reversed because plot and imagesc are opposite, checked this with contour

    title('bar pRF (solid) and SOC pRF (dashed)')
end

plot([res/2 res/2],[0 res],'k')
plot([0 res],[res/2 res/2],'k')
axis off

% save the figure
set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','gaborFilt',...
        ['FigureSupp_BarVSOV']))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','gaborFilt',...
        ['FigureSupp_BarVSOV']))

%%

el_color = jet(length(electrodes));
figure('Position',[0 0 100 200]),hold on
for el = 1:length(electrodes)
    plot(1,el,'o','Color',el_color(el,:),'MarkerSize',10)
end
set(gca,'YTick',1:length(electrodes))
ylim([0 max(length(electrodes))+1])

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300','-painters',fullfile(dataDir,'derivatives','gaborFilt',...
        ['FigureSupp_BarVSOV_colormap']))
print('-dpng','-r300','-painters',fullfile(dataDir,'derivatives','gaborFilt',...
        ['FigureSupp_BarVSOV_colormap']))

