clear all

addpath('~/Documents/git/ecogBasicCode/')
addpath(['/Volumes/DoraBigDrive/data/visual_soc/m-files']);

%% Get electrode labels

clear all
dataRootPath = '/Volumes/DoraBigDrive/data/visual_soc/soc_bids';

subjects = [19,23,24,9];
hemi = {'L','L','R','R'};
hemi_s = {'l','l','r','r'};

v_dirs = [270 0;90 0;90 -60;270 -60;0 0];

for s = 3%1:length(subjects)
    % subject code
    if subjects(s)<10
        subj = ['0' num2str(subjects(s))];
    else
        subj = num2str(subjects(s));
    end
    
    % electrode locations name:
    dataLocName = [dataRootPath '/sub-' subj '/ses-01/ieeg/sub-' subj '_ses-01_acq-corrected_electrodes.tsv'];
    % gifti file name:
    dataGiiName = [dataRootPath '/sub-' subj '/ses-01/anat/sub-' subj '_T1w_pial.' hemi{s} '.surf.gii'];
    % surface labels 
    areas_labels_name = [dataRootPath '/sub-' subj '/ses-01/derivatives/RetinotopyTemplates/rt_sub000/surf/' hemi_s{s} 'h.template_areas.mgz'];
    areas_labels = MRIread(areas_labels_name);
    vert_areas = areas_labels.vol(:);

    angle_labels_name = [dataRootPath '/sub-' subj '/ses-01/derivatives/RetinotopyTemplates/rt_sub000/surf/' hemi_s{s} 'h.template_angle.mgz'];
    angle_labels = MRIread(angle_labels_name);
    vert_angle = angle_labels.vol(:);

    eccen_labels_name = [dataRootPath '/sub-' subj '/ses-01/derivatives/RetinotopyTemplates/rt_sub000/surf/' hemi_s{s} 'h.template_eccen.mgz'];
    eccen_labels = MRIread(eccen_labels_name);
    vert_eccen = eccen_labels.vol(:);
    
    % load electrode locations
    loc_info = readtable(dataLocName,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});
    elecmatrix = [loc_info.x loc_info.y loc_info.z];

    % load gifti:
    g = gifti(dataGiiName);
    
end

%%

% Loop through electrodes and find sphere of closest points on rendering

dist_th = 3; % look 3 mm from electrode
out = [];

elec_label = zeros(size(elecmatrix,1),3); % mode (area), median(angle), median(eccen)

for kk = 1:size(elecmatrix,1)
    a = pdist2(elecmatrix(kk,:),g.vertices,'euclidean');
    out(kk).area = vert_areas(a<dist_th);
    out(kk).angle = vert_angle(a<dist_th);
    out(kk).eccen = vert_eccen(a<dist_th);
end

for kk = 1:size(elecmatrix,1)
    elec_label(kk,1) = mode(out(kk).area);
    elec_label(kk,2) = median(out(kk).angle);
    elec_label(kk,3) = median(out(kk).eccen);
end

save(['./elec_labels/sub-' subj '_elec_label3mm'],'elec_label')