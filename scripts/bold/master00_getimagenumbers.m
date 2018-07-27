%
%
%
% This scripts is used to match fMRI stimulus numbers (from Kay et al.,
% 2013 PLOS Computationl Biology) to image numbers used in ECoG.
%
%
%

% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

%% Load fMRI stimuli

load(fullfile(dataDir,'bold','stimuli.mat'))

%% Write fMRI images as png to compare

% write out the first frame of each stimulus
whframe = 1;  % which index to write out
res = 256;    % downsample to what resolution
for p=1:length(images)
  imwrite(uint8(imresize(double(images{p}(:,:,whframe)),[res res])), ...
          sprintf('thumbnails/thumbnails%03d.png',p));
end

%% Just to check the order of images:

% Display fMRI space stimuli in a figure
i_ind=[70:2:78 79 81:2:83 84 85 86 87:2:99 100 101 102:2:114 115 116,...
    117 118:2:120 122:2:130 131];
figure
for k=1:length(i_ind)
    subplot(6,7,k)
    imagesc(images{i_ind(k)}(:,:,1));
    axis square; colormap gray
    title(['im ' int2str(i_ind(k))])
    ylabel(['new ' int2str(k)])
    set(gca,'XTick',[],'YTick',[])
end

% note: some small differences in some of the space stimuli: from stimulus 7:9 11:12

%%
% Display fMRI contrast stimuli in a figure
i_ind=[139:168];
figure
for k=1:length(i_ind)
    subplot(6,7,k)
    imagesc(images{i_ind(k)}(:,:,1),[0 256]);
    axis square; colormap gray
    title(['im ' int2str(i_ind(k))])
    ylabel(['new ' int2str(38+k)])
    set(gca,'XTick',[],'YTick',[])
end

%%
% Display fMRI separation stimuli in a figure
i_ind=[176:184];
figure
for k=1:length(i_ind)
    subplot(6,7,k)
    imagesc(images{i_ind(k)}(:,:,1),[0 256]);
    axis square; colormap gray
    title(['im ' int2str(i_ind(k))])
    ylabel(['new ' int2str(68+k)])
    set(gca,'XTick',[],'YTick',[])
end

%% Save the conversion: 
% these are the fMRI images that correspond to ECoG images 1:77 (missing 78:86)

imnumbers_ecog=[70:2:78 79 81:2:83 84 85 86 87:2:99 100 101 102:2:114 115 116,...
    117 118:2:120 122:2:130 131 139:168 176:184];
imnumbers_ecog=imnumbers_ecog-69;

save(fullfile(dataDir,'bold','imnumbersbold2ecog'),'imnumbers_ecog')
