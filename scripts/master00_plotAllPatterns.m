% Plot patterns. Original images were developed by Kendrick Kay (PLOS
% Computational Biology, 2013)
%
% Dora Hermes, 2017

clear all
% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath('~/Documents/git/ecogBasicCode/')
addpath(genpath('~/Documents/m-files/knkutils'));

%% load images

load(fullfile(dataDir,'soc_bids','stimuli','task-soc_stimuli.mat'),'stimuli')
stimuli = double(stimuli);

%% figure of all images 

% plot rows of 4 images
for kk = 1:4:86
    f = makeimagestack(stimuli(:,:,kk:kk+3),0,1,[1 4],4);
    figure('Position',[0 0 750 300]);
    imagesc(f);
    axis image tight;
    caxis([0 254]);
    colormap(gray);
    axis off
    set(gcf,'PaperPositionMode','auto')
    print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt',...
        ['stimuli_nr' int2str(kk) '_' int2str(kk+3)]))
    print('-depsc','-r300',fullfile(dataDir,'soc_bids','derivatives','gaborFilt',...
        ['stimuli_nr' int2str(kk) '_' int2str(kk+3)]))
end
