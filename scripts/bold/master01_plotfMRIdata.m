%
%
%
% This scripts is used to plot fMRI beta values (and modelfit) for some
% voxels.
%
%

% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

%% Load fMRI stimuli

% load fMRI stimulus numbers that were used in ECoG
load(fullfile(dataDir,'bold','imnumbersbold2ecog')); % loads imnumbers_ecog

% load fMRI beta values
load(fullfile(dataDir,'bold','dataset03.mat'));
% roi: roi label numbers for each voxel
% roilabels: names corresponding to roi numbers

% load modelfit to fMRI beta values
a = load(fullfile(dataDir,'bold','dataset03_benchmark.mat'));

%% plot the data for a good voxel 

% The following voxel numbers match well with sub 19, electrode 109: 
%    818
%    811
%    647
%    650
%    651 % v2

% find a voxel to plot, either a random good voxel, or a specific number:
% find one:
% % roiInd = 3; % see roilabels
% % fmri_goodVoxels = find(roi==roiInd & vxsselect(:,3)==1);
% % vox_plot = fmri_goodVoxels(1);
% OR give a number
vox_plot = 647;
vox_name = ['v' int2str(vox_plot) 'closetoSub19El109'];

figure('Position',[100 100 600 200]),hold on

% get betas from voxel
bold_plot = betamn(vox_plot,imnumbers_ecog)';
bolderr_plot = betase(vox_plot,imnumbers_ecog)';

% plot
bar(bold_plot,1)
% plot error
plot([1:77'; 1:77'],[bold_plot+bolderr_plot bold_plot-bolderr_plot]','k','LineWidth',1);

% plot stimulus cutoffs
stim_change = [38.5 46.5 50.5 54.5 58.5 68.5 73.5 78.5 82.5];
for k = 1:length(stim_change)
    plot([stim_change(k) stim_change(k)],[0 6.5],'g')
end
text([19 40 47 52 55 61 70 74 79 84],zeros(10,1)+6,{'space','orie','grat','pl','circ','zcon','sp','zsp','coh','nm'})

ax = axis;
axis([0 87 ax(3:4)]);

% %%%%% plot the benchmark perforance on top: %%%%%
% not present for all voxels
% modelfit=a.modelfit;
% socfit_plot = modelfit(vox_plot,imnumbers_ecog(imnumbers_ecog<101));
% plot(1:68,socfit_plot,'r-','LineWidth',3);

set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',[dataDir './bold/figures/BOLD_data_' vox_name])
print('-depsc','-r300',[dataDir './bold/figures/BOLD_data_' vox_name])

%% Plot a voxel location

v_volume = vxs(vox_plot);
zero_vol = zeros(size(meanvol));
zero_vol(v_volume) = 1;
for k=1:size(zero_vol,3)
    [i,j]=find(zero_vol(:,:,k)==1);
    if ~isempty(i)
        ind_vox=[i j k];
    end
end

figure
imagesc(meanvol(:,:,ind_vox(3)))
colormap gray
hold on
% imagesc(zero_vol(:,:,ind_vox(3)))
plot(ind_vox(2),ind_vox(1),'rs')
set(gca,'XTick',[],'YTick',[])

%%set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',['./figures/meanvolwith_' int2str(v_ind)])
% print('-depsc','-r300',['./figures/meanvolwith_' int2str(v_ind)])


