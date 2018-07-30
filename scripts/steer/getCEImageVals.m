function out_vals = getCEImageVals(imEnergy,pp,res,xx,yy)

% Function - getCEImageVals - gets image values from a contrast image: 
% 
% Input: 
%   imEnergy: contrast energy at every pixel, [image X res*res*orientations]
%       can be reshaped to image resXresXorientations by im = reshape(imEnergy(1,:),res,res,8); 
%   pp: x, y, and sigma for the population receptive field
%
% Output:
%   out_vals contains several statistics for the image:
%   out_vals.
%
% D. Hermes 2018

% Average image energies across orientations:
imEnergyMean = zeros(size(imEnergy,1),res*res);
for kk = 1:size(imEnergy,1)
    thisImage = imEnergy(kk,:);
    thisImage = reshape(thisImage,res*res,8);
    imFilt1 = sum(thisImage,2); % Sum across orientations
    imEnergyMean(kk,:) = imFilt1;
end

% Define function for prf center Gaussian:
gaufun1 = @(pp) vflatten(makegaussian2d(res,pp(1),pp(2),pp(3),...
    pp(3),xx,yy,0,0)/(2*pi*pp(3)^2));

% Get values for prf and several increased sizes:
pp1 = gaufun1(pp); %prf
pp2 = gaufun1(pp.* [1 1 2]); % prf 2x size
pp3 = gaufun1(pp.* [1 1 3]); % prf 3x size
pp4 = gaufun1(pp.* [1 1 4]); % prf 4x size

% Make a surround prf:
prfSur = pp2 - .5*pp1; % (2*prf - 1*prf)
prfSur(prfSur<0)=0; % remove values < 0
prfSur = prfSur*(1./sum(prfSur)); % set surface are to 1

% Initiate outputs: Contrast energy values:
out_vals.prf1CE = zeros(size(imEnergyMean,1),1); 
out_vals.prf2CE = zeros(size(imEnergyMean,1),1);
out_vals.prf3CE = zeros(size(imEnergyMean,1),1);
out_vals.prf4CE = zeros(size(imEnergyMean,1),1);
out_vals.totalCE = zeros(size(imEnergyMean,1),1);
out_vals.prfVarCE = zeros(size(imEnergyMean,1),1);
out_vals.prfSurCE = zeros(size(imEnergyMean,1),1);

for kk = 1:size(imEnergyMean,1)
    thisImage = imEnergyMean(kk,:);
    % 1x prf
    out_vals.prf1CE(kk) = sum(thisImage .* pp1');% Sum across space - % imFilt1Sum = the dot product between simple cell contrast energy and weights
    % 2x prf
    out_vals.prf2CE(kk) = sum(thisImage .* pp2');
    % 3x prf
    out_vals.prf3CE(kk) = sum(thisImage .* pp3');
    % 4x prf
    out_vals.prf4CE(kk) = sum(thisImage .* pp4');
    % total CE:
    out_vals.totalCE(kk) = sum(thisImage);
    % variation in CE in the prf:
    out_vals.prfVarCE(kk) = sum(abs((thisImage-sum(thisImage .* pp1'))) .* pp1');
    
    % prf surround CE:
    out_vals.prfSurCE(kk) = sum(thisImage .* prfSur');
end


