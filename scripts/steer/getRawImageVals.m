function out_vals = getRawImageVals(ims,pp,res,xx,yy)

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

% Define function for prf center Gaussian:
gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
    pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

% Get values for prf and several increased sizes:
pp1 = gaufun1(pp); %prf
pp2 = gaufun1(pp.* [1 1 2]); % prf 2x size
pp3 = gaufun1(pp.* [1 1 3]); % prf 3x size
pp4 = gaufun1(pp.* [1 1 4]); % prf 4x size

% Raw image values:
out_vals.minMaxIM = zeros(size(ims,3),1); 
out_vals.varIM = zeros(size(ims,3),1);
out_vals.CEPrfIm = zeros(size(ims,3),1);
out_vals.varPrfIM = zeros(size(ims,3),1);
out_vals.test = zeros(size(ims,3),1);

for kk = 1:size(ims,3)
    im = ims(:,:,kk);
    out_vals.minMaxIM(kk) = max(im(:)) - min(im(:));
    out_vals.varIM(kk) = var(im(:)/mean(ims(:)));
    
    % Mean pixel contrast energy in prf:
    imCon = (im/mean2(im) - 1).^2; % make it a pixel 'contrast energy' image
    imPrf = imCon(:) .* pp1; % imPrf = reshape(im(:) .* gaufun1(pp1),res,res);
    out_vals.CEPrfIm(kk) = sum(imPrf(:)); % var(imPrf(:));

    % Pixel variance in Prf:
    imCon = im/mean2(im) - 1;
    imPrf = imCon(:) .* pp1; % imPrf = reshape(im(:) .* gaufun1(pp1),res,res);
    a = (imCon(:)-sum(imPrf)).^2 .* pp1;
    out_vals.varPrfIM(kk) = sum(a); 
    
    % Test:
    imCon = im;%im/mean2(im) - 1; % could make it a contrast image
    imPrf = imCon - mean(imCon,1); % this introduces the mean vertical contrast in the entire image    
    imPrf = imPrf + fliplr(imPrf); % makes image symmetrical
    imPrf = reshape(imPrf(:) .* pp4,res,res);
    out_vals.test(kk) = sum(imPrf(:).^2);

    % Brightness in prf:
%     out_vals.test(kk) = max(im(:) .* gaufun1(pp1));
%     pp3 = pp1; 
%     pp3(3) = 1*pp1(3);

end

