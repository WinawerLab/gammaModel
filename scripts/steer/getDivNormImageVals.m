function out_vals = getDivNormImageVals(imEnergy,pp,res,xx,yy)

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

n_use = [.5 1 2 4];
r_max = 1;
c50 = .5;

% Initiate outputs: Contrast energy values:
out_vals.NormR = zeros(size(imEnergy,1),length(n_use)); 
out_vals.NormR_n = n_use;
out_vals.NBF = zeros(size(imEnergy,1),1); 

for kk = 1:size(imEnergy,1)
    thisImage = imEnergy(kk,:);
    
    thisImage = reshape(thisImage,res*res,8);
    imFilt1 = thisImage .* repmat(pp2,1,8);
    imFilt1Sum = sum(imFilt1,1); % Sum across space, 
    % imFilt1Sum = the dot product between simple cell contrast energy and weights

    for ii = 1:length(n_use) % different n's to try
        n = n_use(ii);
        
        % Normalization
        dd = (sum(imFilt1Sum.^n));
        nn = sqrt(sum(imFilt1Sum.^2)).^n;
        out_vals.NormR(kk,ii) = r_max * dd .*...
            ((c50.^n + nn).^-1);
%         D_out(kk,ii) = dd;
%         N_out(kk,ii) = nn;
    end

end

% NBF model:
for kk = 1:size(imEnergy,1)
    thisImage = imEnergy(kk,:);
    
    thisImage = reshape(thisImage,res*res,8);
    imFilt1 = thisImage' * pp1;
    out_vals.NBF(kk) = std(imFilt1);
end


