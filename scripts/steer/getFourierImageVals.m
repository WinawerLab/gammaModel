function out_vals = getFourierImageVals(ims,pp,res,xx,yy)

% Function - getRawImageVals - gets image values from black/white image: 
% 
% Input: 
%   ims: images, [res X res X image]
%   pp: x, y, and sigma for the population receptive field
%
% Output:
%   out_vals contains several statistics for the image:
%   out_vals...
%
% D. Hermes 2018

%% TODO: add repeatability within prf and surround, or within hemisphere?

% Define function for prf center Gaussian:
gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
    pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

% Get values for prf and several increased sizes:
pp1 = gaufun1(pp); %prf
pp2 = gaufun1(pp.* [1 1 2]); % prf 2x size
pp3 = gaufun1(pp.* [1 1 3]); % prf 3x size
pp4 = gaufun1(pp.* [1 1 4]); % prf 4x size

% Fourier image values of repetition:
out_vals.rptIm = zeros(size(ims,3),1); 
out_vals.rptPrf = zeros(size(ims,3),1); 

% Fourier for entire image
for kk = 1:size(ims,3)
    im = ims(:,:,kk);    
    im = im/mean2(im) - 1;
       
    a_auto = fftshift(abs(fft2(im)).^2);
    a = mean(mean(a_auto,4),3);

    % Filter image:
%     filt = (fspecial('gaussian',3,3));
%     a=conv2(single(a),filt,'same');
    
    mid_f = floor((size(a,1)+1)/2)+1;
    
    % peaks in first direction:
    a_plot = mean(a,2);
    [PKS,LOCS,w,p] = findpeaks(a_plot,'MinPeakProminence',max(a_plot)/4,'Annotate','extents','WidthReference','halfheight');
    % total power:
    total_power = sum(a_plot(mid_f:end));
    % max power at peaks:
    peak_power = sum(PKS);
    % average width 1:
    av_width(1) = peak_power./total_power;
    
    % peaks in second direction:
    a_plot = mean(a,1);
    [PKS,LOCS,w,p] = findpeaks(a_plot,'MinPeakProminence',max(a_plot)/4,'Annotate','extents','WidthReference','halfheight');
    % total power:
    total_power = sum(a_plot(mid_f:end));
    % max power at peaks:
    peak_power = sum(PKS);
    % average width 1:
    av_width(2) = peak_power./total_power;

    out_vals.rptIm(kk) = (mean(av_width,2)); 
    
    clear av_width

end


% Fourier for X*Prf
for kk = 1:size(ims,3)
    im = ims(:,:,kk);    
    im = im/mean2(im) - 1;
    imPrf = reshape(im(:) .* pp2,res,res); % used to be pp4
   
    a_auto = fftshift(abs(fft2(imPrf)).^2);
    a = mean(mean(a_auto,4),3);

    % middle frequency
    mid_f = floor((size(a,1)+1)/2)+1;
    
    % peaks in first direction:
    a_plot = mean(a,2);
    [PKS,LOCS,w,p] = findpeaks(a_plot,'MinPeakProminence',max(a_plot)/4,'Annotate','extents','WidthReference','halfheight');
    % total power:
    total_power = sum(a_plot(mid_f:end));
    % max power at peaks:
    peak_power = sum(PKS);
    % average width 1:
    av_width(1) = peak_power./total_power;
    
    % peaks in second direction:
    a_plot = mean(a,1);
    [PKS,LOCS,w,p] = findpeaks(a_plot,'MinPeakProminence',max(a_plot)/4,'Annotate','extents','WidthReference','halfheight');
    % total power:
    total_power = sum(a_plot(mid_f:end));
    % max power at peaks:
    peak_power = sum(PKS);
    % average width 1:
    av_width(2) = peak_power./total_power;
    
    out_vals.rptPrf(kk) = (mean(av_width,2)); 
end


%% this is some code that detects peaks in 2D and could get the number of orientations.
% 
% thresh_scale = 15;                      % parameter: how many thresholding steps 
% thresh_perc = 6;                        % parameter: threshold at which we clip
% nhood_minsize = 10;
% 
% thresh = multithresh(a,thresh_scale);    
% q_image = imquantize(a, thresh);         
% 
% q_image(q_image <= thresh_perc) = 0;     % regions under threshold are thrown away
% q_image(q_image > thresh_perc) = 1;      % ... while all others are preserved
% q_image = imbinarize(q_image);           % binarize for further processing
% B = bwareaopen(q_image, nhood_minsize);  % Filter really small regions
% [L, L_num] = bwlabel(B); % <- result     % Label connected components