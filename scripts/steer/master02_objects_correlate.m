
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
%% get PRF info and load ECoG data from one electrode
%%

elec = 68;

% if ismember(elec,[66 67 69])
%     cpfov = 32; % cycles per image
% elseif ismember(elec,[68 104])
%     cpfov = 16; % cycles per image
% elseif ismember(elec,103)
%     cpfov = 8; % cycles per image
% end
   
% Choose an analysis type:
% analysisType = 'spectra200';
% analysisType = 'spectra';
analysisType = 'spectra500';

% load the fitting results:
dataFitName = [dataRootPath '/sub-' subj '/ses-01/derivatives/ieeg/sub-' subj '_task-objects_allruns_' analysisType '_fitEl' int2str(elec) '.mat'];
load(dataFitName)

resamp_parms = resamp_parms_allrepeat;
% resamp_parms = resamp_parms_nonrepeat;

% Broadband estimate, one value per image:
bb_base = resamp_parms(1,6); % from the baseline is the same for resamp_parms(:,6)
ecog_bb = squeeze(resamp_parms(:,2))-bb_base;
ecog_g = squeeze(resamp_parms(:,3));
ecog_a = squeeze(resamp_parms(:,7));

cod_bb = calccod(resamp_parms_evenrepeat(:,2)-bb_base,resamp_parms_oddrepeat(:,2)-bb_base);
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

% Define a model to calculate the contrast within the prf:
modelfun = @(pp,dd) dd*gaufun1(pp);
% modelfun = @(pp,dd) dd*surroundfun(gaufun2(pp),gaufun1(pp));

figure,
imagesc(reshape(gaufun1(pp),res,res),[0 max(gaufun1(pp))]),hold on
contour(xx_1,yy_1,G_el,contourVal,'k')
title('prf')

%% Plot prf features

im_nr_plot = [62];% 50 54 58];

for kk = im_nr_plot%:size(imsR,3) % plot features of this image
    figure('Position',[0 0 800 400])
    thisImage = imEnergy(kk,:);
    thisImage = reshape(thisImage,res*res,8);
    
    imFilt1 = thisImage .* repmat(gaufun1(pp),1,8);
%     imFilt2 = thisImage .* repmat(surroundfun(gaufun2(pp),gaufun1(pp)),1,8);
  
    ax1 = subplot(3,6,1);
    imagesc(ims_repeat(:,:,kk))  
    title(['Stimulus ' int2str(im_nr_plot)])
    axis off
    subplot(3,6,2); 
    imagesc(reshape(mean(thisImage,2),res,res)),hold on
    contour(xx_1,yy_1,G_el,contourVal,'k')
    title('Contrast energy')
    subplot(3,6,3)
    imagesc(reshape(gaufun1(pp),res,res))
%     subplot(3,6,4)
%     imagesc(reshape(surroundfun(gaufun2(pp),gaufun1(pp)),res,res))

    for orientPlot = 1:8
        subplot(3,8,8+orientPlot)
        imagesc(reshape(imFilt1(:,orientPlot),res,res),[0 max(imFilt1(:))])
        axis square, axis off
    end
%     subplot(3,8,12)
%     title('Image within prf filtered by different orientations')
%     for orientPlot = 1:8
%         subplot(3,8,16+orientPlot)
%         imagesc(reshape(imFilt2(:,orientPlot),res,res),[0 max(imFilt2(:))])
%         axis square, axis off
%     end
%     subplot(3,8,20)
%     title('Image surrounding prf filtered by different orientations')
    
    colormap(ax1,'gray')
    
%     set(gcf,'PaperPositionMode','auto')
%     print('-dpng','-r300',['./figures/imageProperties/sub' subj '_el'  int2str(elec) '_image' int2str(kk)])
%     print('-depsc','-r300',['./figures/imageProperties/sub' subj '_el'  int2str(elec) '_image' int2str(kk)])

end

%%
%% Normalization test for all repeated stimuli

pp1 = pp;
% pp1(1) = 220;
% pp1(2) = 150;
% pp1(3) = 5; % electrode 67: / 68:>40, try larger maybe for gamma

im_nr_plot = [1:size(imEnergy,1)];%[6 8 11 52 59 71];% Cars from cerebral cortex: [11 52];

n_use = [.001 .1 .5 1 2 4];
R_out = zeros(length(im_nr_plot),length(n_use));
r_max = 1;
c50 = .7;
% c50 = .7; % electrode 68, try .01 for gamma

for kk = 1:length(im_nr_plot)
    thisImage = imEnergy(im_nr_plot(kk),:);
    thisImage = reshape(thisImage,res*res,8);
    imFilt1 = thisImage .* repmat(gaufun1(pp1),1,8);
    imFilt1Sum = sum(imFilt1,1); % Sum across space
    % imFilt1Sum = the dot product between simple cell contrast energy and weights
    
    for ii = 1:length(n_use)
        n = n_use(ii);
        
        % Normalization
        c_rms = sqrt(sum(imFilt1Sum.^2));
        R_out(kk,ii) = r_max * (sum(imFilt1Sum.^n)).*...
            ((c50.^n + c_rms.^n).^-1);
    end
end

figure
for ii = 1:length(n_use)
    subplot(2,length(n_use),ii)
    plot(R_out(:,ii),ecog_bb,'k.')
    r = corr(R_out(:,ii),ecog_bb);
%     r = corr(R_out(:,ii),ecog_bb,'type','Spearman');
    title(['r=' num2str(r,3)])
    axis tight
end

for ii = 1:length(n_use)
    subplot(2,length(n_use),length(n_use)+ii)
    plot(R_out(:,ii),ecog_g,'k.')
    r = corr(R_out(:,ii),ecog_g);
%     r = corr(R_out(:,ii),ecog_g,'type','Spearman');
    title(['r=' num2str(r,3)])
    axis tight
end


%%
%%
%% Correlate with contrast energy 

pp1 = pp;
% % el 68:
% pp1(1) = 220;
% pp1(2) = 150;
% % el 67:
% pp1(1) = 130;
% pp1(2) = 100;
% pp1(3) = 5; % electrode 67: / 68:>40, try larger maybe for gamma

im_nr_plot = [1:size(imEnergy,1)]; 

imEnergyMean = zeros(length(im_nr_plot),res*res);
for kk = 1:length(im_nr_plot)
    thisImage = imEnergy(im_nr_plot(kk),:);
    thisImage = reshape(thisImage,res*res,8);
    imFilt1 = sum(thisImage,2); % Sum across orientations
    imEnergyMean(kk,:) = imFilt1;
end


%% Contrast energy withing prf

pp1 = pp;

CE_out = zeros(length(im_nr_plot),1);
for kk = 1:length(im_nr_plot)
    thisImage = imEnergyMean(im_nr_plot(kk),:);
    CE_out(kk) = sum(thisImage .* gaufun1(pp1)');% Sum across space - % imFilt1Sum = the dot product between simple cell contrast energy and weights
end

figure('Position',[0 0 600 300])
subplot(1,2,1),hold on
plot(CE_out,ecog_bb,'k.','MarkerSize',10)
%     r = corr(CE_out(:,ii),ecog_bb);
r = corr(CE_out,ecog_bb);%,'type','Spearman');
[b,~,~,~,stats] = regress(ecog_bb,[CE_out ones(size(CE_out))]);
x = [0:0.1:1.5];
plot(x,b(1)*x+b(2),'k')
title(['r=' num2str(r,3)])
xlim([min(CE_out-.1) max(CE_out+.1)]),ylim([min(ecog_bb)-.01 max(ecog_bb)+.01])
plot(CE_out([11]),ecog_bb([11]),'b.','MarkerSize',20)
plot(CE_out([52]),ecog_bb([52]),'r.','MarkerSize',20)
xlabel('contrast energy in gaufun(pp)'),ylabel('broadband')

subplot(1,2,2),hold on
plot(CE_out,ecog_g,'k.','MarkerSize',10)
%     r = corr(CE_out(:,ii),ecog_g);
r = corr(CE_out,ecog_g);%,'type','Spearman');
% [b,~,~,~,stats] = regress(ecog_g,[CE_out ones(size(CE_out))]);
[b] = polyfit(CE_out,ecog_g,1);
x = [0:0.1:1.5];
plot(x,b(1)*x+b(2),'k')
title(['r=' num2str(r,3)])
% % xlim([min(CE_out-.1) max(CE_out+.1)]),ylim([min(ecog_g)-.001 max(ecog_g)+.001])
plot(CE_out([11]),ecog_g([11]),'b.','MarkerSize',20)
plot(CE_out([52]),ecog_g([52]),'r.','MarkerSize',20)
xlabel('contrast energy in gaufun(pp)'),ylabel('gamma')

%% Correlate contrast energy with ECoG for every pixel

corr_bb = zeros(1,size(imEnergyMean,2));
corr_g = zeros(1,size(imEnergyMean,2));
for kk = 1:size(imEnergyMean,2)
    if mod(kk,1000)==0
        disp(['pixel ' int2str(kk) ' of '  int2str(size(imEnergyMean,2))])
    end
    corr_bb(kk) = corr(imEnergyMean(:,kk),ecog_bb);
    corr_g(kk) = corr(imEnergyMean(:,kk),ecog_g);
end

figure('Position',[0 0 600 300])
subplot(1,2,1)
imagesc(reshape(corr_bb.^2,res,res),[0 .4]),hold on
contour(xx_1,yy_1, G_el,contourVal,'k')
axis square, axis off
colorbar
title('r^2: CE and bb')

subplot(1,2,2)
imagesc(reshape(corr_g.^2,res,res),[0 .4]),hold on
contour(xx_1,yy_1, G_el,contourVal,'k')
axis square, axis off
colorbar
title('r^2: CE and gamma')

set(gcf,'PaperPositionMode','auto')
% print('-depsc','-r300',['./figures/corr_CE_ECoG/sub-' subj '_' analysisType '_el' int2str(elec) '_corrEcogCe'])
% print('-dpng','-r300',['./figures/corr_CE_ECoG/sub-' subj '_' analysisType '_el' int2str(elec) '_corrEcogCe'])

%%
%% Correlate variance across orientations ECoG for every pixel

imEnergyVar = std(reshape(imEnergy,size(imEnergy,1),[],8),[],3);

corr_bb = zeros(1,size(imEnergyVar,2));
corr_g = zeros(1,size(imEnergyVar,2));
for kk = 1:size(imEnergyVar,2)
    if mod(kk,1000)==0
        disp(['pixel ' int2str(kk) ' of '  int2str(size(imEnergyVar,2))])
    end
    corr_bb(kk) = corr(imEnergyVar(:,kk),ecog_bb);
    corr_g(kk) = corr(imEnergyVar(:,kk),ecog_g);
end

figure('Position',[0 0 600 300])
subplot(1,2,1)
imagesc(reshape(corr_bb.^2,res,res),[0 .2]),hold on
contour(xx_1,yy_1, G_el,contourVal,'k')
axis square, axis off
colorbar
title('r^2: NBF and bb')

subplot(1,2,2)
imagesc(reshape(corr_g.^2,res,res),[0 .2]),hold on
contour(xx_1,yy_1, G_el,contourVal,'k')
axis square, axis off
colorbar
title('r^2: NBF and gamma')

set(gcf,'PaperPositionMode','auto')
print('-depsc','-r300',['./figures/corr_CE_ECoG/sub-' subj '_' analysisType '_el' int2str(elec) '_corrEcogNBF'])
print('-dpng','-r300',['./figures/corr_CE_ECoG/sub-' subj '_' analysisType '_el' int2str(elec) '_corrEcogNBF'])

%% Normalization test for all repeated stimuli - with Gating Term

pp1 = pp;
% pp1(1) = 220;
% pp1(2) = 150;
% pp1(3) = 5; % electrode 67/68, try 20 for gamma

im_nr_plot = [1:size(imEnergy,1)];%[6 8 11 52 59 71];% Cars from cerebral cortex: [11 52];

n_use = [.1 .5 1 2 4];
R_out = zeros(length(im_nr_plot),length(n_use));
r_max = 1;
c50 = .7;
% c50 = .2; % electrode 68, try .01 for gamma

prf_center = repmat(gaufun1(pp1),1,8);
prf_surround1 = repmat(gaufun1(pp1 + [1.7*pp1(3) 0 0]),1,8);
prf_surround2 = repmat(gaufun1(pp1 + [1.7*pp1(3) 1.7*pp1(3) 0]),1,8);
prf_surround3 = repmat(gaufun1(pp1 + [0 1.7*pp1(3) 0]),1,8);
prf_surround4 = repmat(gaufun1(pp1 + [-1.7*pp1(3) 1.7*pp1(3) 0]),1,8);
prf_surround5 = repmat(gaufun1(pp1 + [-1.7*pp1(3) 0 0]),1,8);
prf_surround6 = repmat(gaufun1(pp1 + [-1.7*pp1(3) -1.7*pp1(3) 0]),1,8);
prf_surround7 = repmat(gaufun1(pp1 + [0 -1.7*pp1(3) 0]),1,8);
prf_surround8 = repmat(gaufun1(pp1 + [1.7*pp1(3) -1.7*pp1(3) 0]),1,8);
prf_surround = cat(3,prf_surround1,prf_surround2,prf_surround3,prf_surround4,prf_surround5,prf_surround6,prf_surround7,prf_surround8);

gatingTerm = zeros(length(im_nr_plot),1);

for kk = 1:length(im_nr_plot)
    thisImage = imEnergy(im_nr_plot(kk),:);
    thisImage = reshape(thisImage,res*res,8);
    
    % Center pool
    imFilt1 = thisImage .* prf_center; % imFilt1Sum = the dot product between simple cell contrast energy and weights
    imFilt1Sum = sum(imFilt1,1); % Sum across space

    % Surround pool
    imFilt2Sum = zeros(8,size(prf_surround,3));
    for ss = 1:size(prf_surround,3)
        imFilt2 = thisImage .* prf_surround(:,:,ss);
        imFilt2Sum(:,ss) = sum(imFilt2,1); % Sum across space
    end
    
    % Gating term - image homogeneity
    % this term should be higher for more homogeneous patterns
    center_surround_homo = zeros(size(prf_surround,3),1);
    for ss = 1:size(prf_surround,3)
        center_surround_homo(ss) = sum(imFilt1Sum.*imFilt2Sum(:,ss)') ./ (norm(imFilt1Sum)*norm(imFilt2Sum(:,ss)')); % removes variation in contrast
    end
    gatingTerm(kk) = mean(center_surround_homo);    
    
    imFilt2Sum = mean(imFilt2Sum,2);
    
    for ii = 1:length(n_use)
        n = n_use(ii);
        
        % Normalization
        c_rms = sqrt(sum(imFilt2Sum.^2));
        R_out(kk,ii) = r_max * (sum(imFilt1Sum.^n)).*...
            ((c50.^n + c_rms.^n).^-1);
    end
end

figure
for ii = 1:length(n_use)
    subplot(2,length(n_use),ii)
    plot(R_out(:,ii),ecog_bb,'k.')
    r = corr(R_out(:,ii),ecog_bb);
%     r = corr(R_out(:,ii),ecog_bb,'type','Spearman');
    title(['r=' num2str(r,3)])
    axis tight
end

for ii = 1:length(n_use)
    subplot(2,length(n_use),length(n_use)+ii)
    plot(R_out(:,ii),ecog_g,'k.')
    r = corr(R_out(:,ii),ecog_g);
%     r = corr(R_out(:,ii),ecog_g,'type','Spearman');
    title(['r=' num2str(r,3)])
    axis tight
end

%% now, what is the difference between these images:
% 2 groups:
% 1) high contrast and low gamma % is this image heterogeneous?
% 2) high contrast and high gamma
% Question:
% higher gating term --> more homogeneous --> more gamma?
% no: opposite: corr(gatingTerm,ecog_g) is negative (?)

gatingTerm = gatingTerm;

figure,hold on
% these are some images with high contrast energy, but low gamma:
test_ims1 = find(CE_out>.8 & ecog_g<.02);
bar(1,mean(gatingTerm(test_ims1)),'c')
% plot(1,gatingTerm(test_ims1),'k.') 
errorbar(1,mean(gatingTerm(test_ims1)),std(gatingTerm(test_ims1))./sqrt(length(test_ims1)),'k')

% these are some images with high contrast energy, and high gamma:
test_ims2 = find(CE_out>.8 & ecog_g>.04);
bar(2,mean(gatingTerm(test_ims2)),'c')
% plot(2,gatingTerm(test_ims2),'k.') 
errorbar(2,mean(gatingTerm(test_ims2)),std(gatingTerm(test_ims2))./sqrt(length(test_ims2)),'k')

xlim([0 3])

[~,p,~,~] = ttest2(gatingTerm(test_ims1),gatingTerm(test_ims2))

%%
%% FFT for periodicity in prf and surround 

pp1 = pp;
pp1(3) = 4*pp1(3);

% Define outputs:
av_width = zeros(size(ims_repeat,3),2);
im_test = zeros(size(ims_repeat,3),1);
% nr_orientations = zeros(size(ims_repeat,3),1); % TODO, only makes sense if repetitive...

for kk = 1:size(ims_repeat,3)
    im = double(ims_repeat(:,:,kk));
    
    im_test(kk) = mean(im(:));
    
    im = im/mean2(im) - 1;
    imPrf = reshape(im(:) .* gaufun1(pp1),res,res);
   
    a_auto = fftshift(abs(fft2(imPrf)).^2);
    a = mean(mean(a_auto,4),3);
    
    % Filter in Fourier or not, may be better for estimating the number of
    % similarity-axes/orientations
    filt = (fspecial('gaussian',3,3));
    a=conv2(single(a),filt,'same');

    mid_f = floor((size(a,1)+1)/2)+1;
    
    % peaks in first direction:
    a_plot = mean(a,2);
    [PKS,LOCS,w,p] = findpeaks(a_plot,'MinPeakProminence',max(a_plot)/4,'Annotate','extents','WidthReference','halfheight');
    % total power:
    total_power = sum(a_plot(mid_f:end));
    % max power at peaks:
    peak_power = sum(PKS);
    % average width 1:
    av_width(kk,1) = peak_power./total_power;
    
    % peaks in second direction:
    a_plot = mean(a,1);
    [PKS,LOCS,w,p] = findpeaks(a_plot,'MinPeakProminence',max(a_plot)/4,'Annotate','extents','WidthReference','halfheight');
    % total power:
    total_power = sum(a_plot(mid_f:end));
    % max power at peaks:
    peak_power = sum(PKS);
    % average width 1:
    av_width(kk,2) = peak_power./total_power;
    
end

level_period = mean(av_width,2);

figure,plot(level_period)

% Contrast energy withing prf
pp1 = pp;
im_nr_plot = [1:size(imEnergy,1)]; 
imEnergyMean = zeros(length(im_nr_plot),res*res);
for kk = 1:length(im_nr_plot)
    thisImage = imEnergy(im_nr_plot(kk),:);
    thisImage = reshape(thisImage,res*res,8);
    imFilt1 = sum(thisImage,2); % Sum across orientations
    imEnergyMean(kk,:) = imFilt1;
end

CE_out = zeros(length(im_nr_plot),1);
for kk = 1:length(im_nr_plot)
    thisImage = imEnergyMean(im_nr_plot(kk),:);
    CE_out(kk) = sum(thisImage .* gaufun1(pp1)');% Sum across space - % imFilt1Sum = the dot product between simple cell contrast energy and weights
end

%%
% im_noisy = level_period<median(level_period);
% im_noisy1 = level_period<quantile(level_period,.13);
% im_noisy2 = level_period<quantile(level_period,1-.13);

[B,BINT,R,RINT,STATS]  = regress(ecog_g,[CE_out ones(size(CE_out))]);
STATS(1)
[B,BINT,R,RINT,STATS]  = regress(ecog_g,[CE_out im_noisy1 im_noisy2 ones(size(CE_out))]);
STATS(1)


%% test crazy stuff that works (??)

pp1 = pp;
% pp1(1) = 175; % sanity check
% pp1(2) = 175; % sanity check
pp1(3) = 4*pp1(3);

im_out.min_max = zeros(size(ims_repeat,3),1);
im_out.min_maxPrf = zeros(size(ims_repeat,3),1);
im_out.var = zeros(size(ims_repeat,3),1);
im_out.varPrf = zeros(size(ims_repeat,3),1);
for kk = 1:size(ims_repeat,3)
    im = ims_repeat(:,:,kk);
    im_out.min_max(kk) = max(im(:)) - min(im(:));
    im_out.var(kk) = var(im(:));
    
% for some reason, this correlates with gamma, subtracting the mean is important and the calculating the variance:
    imCon = im/mean2(im) - 1; % make it a contrast image
    imPrf = imCon - mean(imCon,1); % this introduces the mean vertical contrast in the entire image    
    imPrf = imPrf + fliplr(imPrf); % makes image symmetrical
    
    imPrf = reshape(imPrf(:) .* gaufun1(pp1),res,res);
%     im_out.varPrf(kk) = var(imPrf(:));
    im_out.varPrf(kk) = sum(imPrf(:).^2);

% this correlates less well with gamma:
%     imPrf = im(:) .* gaufun1(pp1);
%     imPrf(gaufun1(pp1)<.3*max(gaufun1(pp1))) = [];
%     imPrf = (im/mean(imPrf) - 1) - mean(imPrf);
%     im_out.varPrf(kk) = var(imPrf(:));
end

corr(im_out.min_max,ecog_g).^2
corr(im_out.var,ecog_g).^2
corr(im_out.varPrf,ecog_g).^2

%% test crazy stuff that works (??) on ims_nonrepeat

pp1 = pp;
pp1(3) = 4*pp1(3);

im_out.min_max = zeros(size(ims_nonrepeat,3),1);
im_out.min_maxPrf = zeros(size(ims_nonrepeat,3),1);
im_out.var = zeros(size(ims_nonrepeat,3),1);
im_out.varPrf = zeros(size(ims_nonrepeat,3),1);
for kk = 1:size(ims_nonrepeat,3)
    im = ims_nonrepeat(:,:,kk);
    im_out.min_max(kk) = max(im(:)) - min(im(:));
    im_out.var(kk) = var(im(:));
    
% for some reason, this correlates with gamma, subtracting the mean is important and the calculating the variance:
    imCon = im/mean2(im) - 1; % make it a contrast image
    imPrf = imCon - mean(imCon,1); % this introduces the mean vertical contrast in the entire image    
    imPrf = imPrf + fliplr(imPrf); % makes image symmetrical
    
    imPrf = reshape(imPrf(:) .* gaufun1(pp1),res,res);
%     im_out.varPrf(kk) = var(imPrf(:));
%     im_out.varPrf(kk) = sum(imPrf(:).^2);
    im_out.varPrf(kk) = sum(abs(imPrf(:)));

% this correlates less well with gamma:
%     imPrf = im(:) .* gaufun1(pp1);
%     imPrf(gaufun1(pp1)<.3*max(gaufun1(pp1))) = [];
%     imPrf = (im/mean(imPrf) - 1) - mean(imPrf);
%     im_out.varPrf(kk) = var(imPrf(:));
end

corr(im_out.min_max,ecog_g)
corr(im_out.var,ecog_g)
corr(im_out.varPrf,ecog_g)
% corr(CE_out,ecog_g)


