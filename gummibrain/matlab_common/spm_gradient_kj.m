function Q = spm_gradient_kj(P,Q,s,graddim)
% 3 dimensional convolution of an image
% FORMAT spm_smooth(P,Q,S,dtype)
% P     - image to be smoothed (or 3D array)
% Q     - filename for smoothed image (or 3D array)
% S     - [sx sy sz] Gaussian filter width {FWHM} in mm (or edges)
% dtype - datatype [default: 0 == same datatype as P]
%____________________________________________________________________________
%
% spm_smooth is used to smooth or convolve images in a file (maybe).
%
% The sum of kernel coeficients are set to unity.  Boundary
% conditions assume data does not exist outside the image in z (i.e.
% the kernel is truncated in z at the boundaries of the image space). S
% can be a vector of 3 FWHM values that specifiy an anisotropic
% smoothing.  If S is a scalar isotropic smoothing is implemented.
%
% If Q is not a string, it is used as the destination of the smoothed
% image.  It must already be defined with the same number of elements
% as the image.
%
%_______________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% John Ashburner & Tom Nichols
% $Id: spm_smooth.m 4172 2011-01-26 12:13:29Z guillaume $


%-----------------------------------------------------------------------
if numel(s) == 1, s = [s s s]; end
if nargin < 4, graddim = 1; end;

if ischar(P), P = spm_vol(P); end;

if isstruct(P),
    Q0 = Q;
    for i= 1:numel(P),
        Q(i) = smooth1(P(i),Q0,s,graddim);
    end
else
    Q = smooth1(P,Q,s,graddim);
end
%_______________________________________________________________________

%_______________________________________________________________________
function Q = smooth1(P,Q,s,graddim)
if isstruct(P),
    VOX = sqrt(sum(P.mat(1:3,1:3).^2));
else
    VOX = [1 1 1];
end;

if ischar(Q) && isstruct(P),
    [pth,nam,ext,num] = spm_fileparts(Q);
    q         = fullfile(pth,[nam,ext]);
    Q         = P;
    Q.fname   = q;
    if ~isempty(num),
        Q.n       = str2num(num);
    end;
    if ~isfield(Q,'descrip'), Q.descrip = sprintf('SPM compatible'); end;
    Q.descrip = sprintf('%s - gradient(dim %d), size(%g,%g,%g)',Q.descrip, graddim,s);

end

% compute parameters for spm_conv_vol
%-----------------------------------------------------------------------
s  = s./VOX;                        % voxel anisotropy
s1 = s/sqrt(8*log(2));              % FWHM -> Gaussian parameter

d = 10;
x  = round(6*s1(1)); x = -x:1/d:x; x = spm_smoothkern(s(1),x,1); x  = x/sum(x);
y  = round(6*s1(2)); y = -y:1/d:y; y = spm_smoothkern(s(2),y,1); y  = y/sum(y);
z  = round(6*s1(3)); z = -z:1/d:z; z = spm_smoothkern(s(3),z,1); z  = z/sum(z);

[xi yi zi] = meshgrid(1:numel(x),1:numel(y),1:numel(z));
S = x(xi).*y(yi).*z(zi);

[Sx Sy Sz] = gradient(S);

Sx = Sx(1:d:end,1:d:end,1:d:end);
Sy = Sy(1:d:end,1:d:end,1:d:end);
Sz = Sz(1:d:end,1:d:end,1:d:end);

Sx = Sx/sqrt(sum(Sx(:).^2));
Sy = Sy/sqrt(sum(Sy(:).^2));
Sz = Sz/sqrt(sum(Sz(:).^2));

Sxyz = {Sx,Sy,Sz};
G = cell(1,numel(graddim));
for i = 1:numel(graddim)
    G{i} = conv3d(P,Sxyz{graddim(i)});
end
G = cat(4,G{:});

if isstruct(Q)
    Q = spm_create_vol(Q);
    spm_write_vol(Q,G);
else
    Q = G;
end

% 
% x = gradient(x);
% y = gradient(y);
% z = gradient(z);
% 
% x = x(1:d:end);
% y = y(1:d:end);
% z = z(1:d:end);
% 
% i  = (length(x) - 1)/2;
% j  = (length(y) - 1)/2;
% k  = (length(z) - 1)/2;
% 
% 
% if isstruct(Q), Q = spm_create_vol(Q); end;
% spm_conv_vol(P,Q,x,y,z,-[i,j,k]);
% 
% 
