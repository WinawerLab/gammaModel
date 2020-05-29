
% This script renders the cortices with or without visual maps with
% electrodes.
%
% dhermes 2018 UMC Utrecht

% set paths:
gammaModelCodePath;
dataDir = gammaModelDataPath;

% add other toolboxes:
addpath('/Users/m206305/Documents/git/ecogBasicCode/render/')
addpath(genpath('/Users/m206305/Documents/git/knkutils'));

%% Write gifti in the space of the original MRI:

subjects = {'19','23','24','9','1001'};
hemi = {'L','L','R','R','R'};
hemi_2 = {'lh','lh','rh','rh','rh'};

% which subject to process
s = 5;

subj = subjects{s};

%%% load the gifti file created from the freesurfer pial by mris_convert lh.pial lh.pial.gii
g = gifti(fullfile(dataDir,'derivatives','freesurfer',['sub-' subj],'surf',...
    [hemi_2{s} '.pial.gii']));

% covert to original space:
mri_orig = fullfile(dataDir,'derivatives','freesurfer',['sub-' subj],'mri','orig.mgz');
orig = MRIread(mri_orig);
Torig = orig.tkrvox2ras;
Norig = orig.vox2ras;
freeSurfer2T1 = Norig*inv(Torig);

% convert vertices to original space
vert_mat = double(([g.vertices ones(size(g.vertices,1),1)])');
vert_mat = freeSurfer2T1*vert_mat;
vert_mat(4,:) = [];
vert_mat = vert_mat';
g.vertices = vert_mat; clear vert_mat

% save as a gifti
gifti_name = fullfile(dataDir,'derivatives','render',['sub-' subj],...
    ['sub-' subj '_T1w_pial.' hemi{s} '.surf.gii']);

save(g,gifti_name,'Base64Binary')


%% Render plain with used electrodes
clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = {'19','23','24','9','1001'};
hemi = {'L','L','R','R','R'};
els_gammaModel = {[107 108 109 115 120 121],'',[45 46],'',[49 50 52 57 58 59 60]};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

for s = 5%1:length(subjects)
    % subject code    
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','render',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);

    % electrode locations name:
    dataLocName = fullfile(dataRootPath,['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_acq-corrected_electrodes.tsv']);
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];

    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);

        figure
        ecog_RenderGifti(g) % render
%         ecog_Label(els,30,12) % add electrode positions
        el_add(els(els_gammaModel{s},:),'k',30) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        
%         set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',fullfile(dataDir,'soc_bids','derivatives','render',['sub-' subj],...
%             ['subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))]))

        % close all
    end
end


%% Render plain with all electrodes

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = {'19','23','24','9','1001'};
hemi = {'L','L','R','R','R'};

v_dirs = [270 0;90 0;90 -60;270 -60;90 0];
v_dirs = [270 0;90 0;90 -60;270 -60;-60 -15];

for s = 5%1:length(subjects)
    % subject code  
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','render',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);

    % electrode locations name:
    dataLocName = fullfile(dataRootPath,['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_acq-corrected_electrodes.tsv']);
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];

    
    % figure with rendering for different viewing angles
    v_d = v_dirs(s,:);

    % make sure electrodes pop out
    a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
    els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);

    figure
    ecog_RenderGifti(g) % render
    ecog_Label(els,30,12) % add electrode positions
    ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
    set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['../figures/render/subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
%         close all
end

%% Render Wang/Kastner with all electrodes

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = {'19','23','24','9','1001'};
hemi = {'L','L','R','R','R'};
hemi_s = {'l','l','r','r'};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','render',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);

    % electrode locations name:
    dataLocName = fullfile(dataRootPath,['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_acq-corrected_electrodes.tsv']);
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % surface labels 
    surface_labels_name = fullfile(dataRootPath,'derivatives','freesurfer',['sub-' subj],'surf',...
        [hemi_s{s} 'h.wang15_mplbl.mgz']);
    surface_labels = MRIread(surface_labels_name);
    vert_label = surface_labels.vol(:);

    % cmap = 'lines';
    cmap = lines(max(vert_label));
    
    % figure with rendering for different viewing angles
    for k = 1%:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        % make sure electrodes pop out
%         a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
%         els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);
        els = elecmatrix;

        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,Wang_ROI_Names)
        
%         ecog_RenderGifti(g) % render
        ecog_Label(els,30,12) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['./figures/render/Wang_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
%         close all
    end
end


%% Render Benson Areas with electrodes

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids_sourcedir';

subjects = {'19','23','24','9','1001'};
hemi = {'L','L','R','R','R'};
hemi_s = {'l','l','r','r','r'};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0;0 -90;22 -20;-22 -20];

Benson_Area_Names = {'V1','V2','V3','hV4','V01','V02','V3a','V3b','LO1','LO2','TO1','T02'};

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','render',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);

    % electrode locations name:
    dataLocName = fullfile(dataRootPath,['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_acq-corrected_electrodes.tsv']);
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % surface labels 
    surface_labels_name = fullfile(dataRootPath,'derivatives','freesurfer',['sub-' subj],'surf',...
        [hemi_s{s} 'h.benson14_varea.mgz']);
    surface_labels = MRIread(surface_labels_name);
    vert_label = surface_labels.vol(:);

    % cmap = 'lines';
%     cmap = lines(max(vert_label));
    cmap = make_BensonColormap;
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);

%         els = elecmatrix;
        
        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,Benson_Area_Names)
        
%         ecog_RenderGifti(g) % render
        ecog_Label(els,30,12) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
        figureName = fullfile(dataRootPath,'derivatives','render',['sub-' subj],...
            ['sub-' subj '_BensonAreas_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2)) 'allLabels']);

        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',figureName)
        close all
    end
end
   


%%
%% Render Benson Eccentricity with electrodes

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = {'19','23','24','9','1001'};
hemi = {'L','L','R','R','R'};
hemi_s = {'l','l','r','r'};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','render',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);

    % electrode locations name:
    dataLocName = fullfile(dataRootPath,['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_acq-corrected_electrodes.tsv']);
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % surface labels 
    surface_labels_name = fullfile(dataRootPath,'derivatives','freesurfer',['sub-' subj],'surf',...
        [hemi_s{s} 'h.benson14_eccen.mgz']);
    surface_labels = MRIread(surface_labels_name);
    vert_label = surface_labels.vol(:);
    
    % colormap:
    [cm,cm_angle,cm_eccen] = make_BensonColormap();
    cmap = cm_eccen; % or try a standard map like hsv(ceil(max(vert_label)));
    Benson_Eccen_Names = [1:ceil(max(vert_label))];
    
    % multiply vert_label times 4 to use range from 1-90 to 1-360 to use
    % the entire colormap
    vert_label = vert_label*4;
    
    % figure with rendering for different viewing angles
    for k = 1%:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);

        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,Benson_Eccen_Names)
        
%         ecog_RenderGifti(g) % render
        ecog_Label(els,30,12) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
%         set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['./figures/render/BensonEccen_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
%         close all
    end
end


%% Render Benson Angle with electrodes

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = {'19','23','24','9','1001'};
hemi = {'L','L','R','R','R'};
hemi_s = {'l','l','r','r'};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

for s = 1%1:length(subjects)
    % subject code
    subj = subjects{s};
    
    % gifti file name:
    dataGiiName = fullfile(dataRootPath,'derivatives','render',['sub-' subj],...
        ['sub-' subj '_T1w_pial.' hemi{s} '.surf.gii']);
    % load gifti:
    g = gifti(dataGiiName);

    % electrode locations name:
    dataLocName = fullfile(dataRootPath,['sub-' subj],'ses-01','ieeg',...
        ['sub-' subj '_ses-01_acq-corrected_electrodes.tsv']);
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];
    
    % surface labels 
    surface_labels_name = fullfile(dataRootPath,'derivatives','freesurfer',['sub-' subj],'surf',...
        [hemi_s{s} 'h.benson14_angle.mgz']);
    surface_labels = MRIread(surface_labels_name);
    vert_label = surface_labels.vol(:);
    
    [cm,cm_angle,cm_eccen] = make_BensonColormap();
    cmap = cm_angle;%hsv(ceil(max(vert_label)));
    Benson_Angle_Names = [1:ceil(max(vert_label))];
    
    % figure with rendering for different viewing angles
    for k = 1:size(v_dirs,1) % loop across viewing angles
        v_d = v_dirs(k,:);
        
        % make sure electrodes pop out
        a_offset = .1*max(abs(elecmatrix(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = elecmatrix+repmat(a_offset,size(elecmatrix,1),1);

        figure
        ecog_RenderGiftiLabels(g,vert_label,cmap,Benson_Angle_Names)
        
%         ecog_RenderGifti(g) % render
        ecog_Label(els,30,12) % add electrode positions
        ecog_ViewLight(v_d(1),v_d(2)) % change viewing angle   
%         set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',['./figures/render/BensonAngle_subj_' subj '_v' int2str(v_d(1)) '_' int2str(v_d(2))])
%         close all
    end
end


