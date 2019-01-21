function cm = make_BensonColormap()

% This colormap was made using linspecer:
% Beautiful and distinguishable line colors + colormap
% version 1.4.0.0 (8.25 KB) by Jonathan C. Lansey
% and then made darker/lighter for a cortical surface
% and then shuffled to make sure adjacent maps differ in color

% Benson_Area_Names = {'V1','V2','V3','hV4','V01','V02','V3a','V3b','LO1','LO2','TO1','T02'};

cm = [...
    0.9438    0.3910    0.2668
    0.1974    0.5129    0.7403
    0.5978    0.8408    0.6445
    0.9955    0.8227    0.4828
    0.3686    0.3098    0.6353
    0.8417    0.9409    0.6096
    0.6196    0.0039    0.2588
    0.8    0.8    0.4
    0.3539    0.7295    0.6562
    0.9877    0.6154    0.3391
    0.8206    0.2239    0.3094
    0.9    0.6    1];



% % Benson ELife 2018 colors areas:
% cm = [...
%     1,0,0;...%V1
%     0,1,0;...%V2
%     0,0,1;...%V3
%     0,1,1;...%hV4
%     0,0.5,1;...%V01
%     0,1,0.5;...%VO2
%     1,1,0;...'V3a'
%     0.5,1,0;...'V3b'
%     1,0,1;...'LO1'
%     0.5,0,1;...'LO2'
%     1,0,0.5;...'TO1'
%     0.5,0,0.5];%'TO2'
    

% Benson 2018 ELife colors angle
% cmap_polar_angle is a colormap for plotting the pRF polar angle of a vertex.
% Values passed to cmap_polar_angle should be scaled such that (-180,180 deg) -> (0,1).
cmap_polar_angle = [...
    0.5,0,0;...
    1,0,1;...
    0,0,0.5;...
    0,1,1;...
    0,0.5,0;...
    1,1,0;...
    0.5,0,0];
cm_angle = zeros(180,3);
cm_angle(1:)
-180°,-135°,-90°,-45°,0°,45°,90°,135°,180° 


% Benson 2018 ELife colors eccentricity 

cmap_eccentricity = blend_cmap(
        'eccentricity',
        [(0,       (  0,  0,  0)),
         (1.25/90, (  0,  0,0.5)),
         (2.5/90,  (  1,  0,  1)),
         (5.0/90,  (0.5,  0,  0)),
         (10.0/90, (  1,  1,  0)),
         (20.0/90, (  0,0.5,  0)),
         (40.0/90, (  0,  1,  1)),
         (1,       (  1,  1,  1))])
    cmap_eccentricity.__doc__ = '''
    cmap_eccentricity is a colormap for plotting the pRF eccentricity of a vertex.
    Values passed to cmap_eccentricity should be scaled such that (0,90 deg) -> (0,1).
    '''
    cmap_log_eccentricity = blend_cmap(
        'log_eccentricity',
        [(0,     (  0,  0,  0)),
         (1.0/7, (  0,  0,0.5)),
         (2.0/7, (  1,  0,  1)),
         (3.0/7, (0.5,  0,  0)),
         (4.0/7, (  1,  1,  0)),
         (5.0/7, (  0,0.5,  0)),
         (6.0/7, (  0,  1,  1)),
         (1,     (  1,  1,  1))])
    cmap_log_eccentricity.__doc__ = '''
    cmap_log_eccentricity is a colormap for plotting the log of eccentricity.
    Values passed to cmap_log_cmag should be scaled however desired, but note that the colormap
    itself runs linearly from 0 to 1, so eccentricity data should be log-transformed
    then scaled before being passed.
             