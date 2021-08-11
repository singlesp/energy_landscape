function img = Colour_Wheel(wheel, start, cycles, direc)
%%
% img = Colour_Wheel(wheel, start, cycles, direc)
%
% Generates the polar angle RYGB colour wheel for Freesurfer.
%   
%   wheel :     'Polar' or 'Eccen'
%
%   start :     Polar angle or position where cycling begins (0=3 o'clock, 90=12 o'clock)
%
%   cycles :    Number of cycles for the colour wheel. This is the same as the 
%               "angle cycles" defined in the Freesurfer options. You must use 
%               an integer here. Freesurfer does funky stuff when using fractions 
%               for this value, but this script does not support it.
%
%   direc :     Direction of cycling: 
%                   0 = Clockwise rotating wedge / Expanding ring
%                   1 = Anticlockwise rotating wedge / Contracting ring
%
% The output of this function is a bitmap which can be displayed or saved.
%   (Written by Sam Schwarzkopf, 2nd June 2010
% 

radius = 200;    % radius of colour wheel
steps = 360 / 4 / cycles;   % steps from one colour to the next

% colour peaks
R = [1 0 0];
Y = [1 1 0];
G = [0 1 0];
B = [0 0 1];

% create the colour map
cmap = [];
for c = 1:cycles
    cmap = [cmap; linspace(R(1), Y(1), steps)', linspace(R(2), Y(2), steps)', linspace(R(3), Y(3), steps)'];
    cmap = [cmap; linspace(Y(1), G(1), steps)', linspace(Y(2), G(2), steps)', linspace(Y(3), G(3), steps)'];
    cmap = [cmap; linspace(G(1), B(1), steps)', linspace(G(2), B(2), steps)', linspace(G(3), B(3), steps)'];
    cmap = [cmap; linspace(B(1), R(1), steps)', linspace(B(2), R(2), steps)', linspace(B(3), R(3), steps)'];
end

% image pixels
[x y] = meshgrid(-radius:radius, -radius:radius);
if ~direc
    x = flipud(x); y = flipud(y);   % flipping upside down because of the way matlab matrices work
end
[t r] = cart2pol(x,y);
t = t / pi * 180;   % convert to degrees
t = t + 90 + start;  % calibrate angles to stimulus starting position
t = mod(ceil(t),360) + 1;   % ensure between 1-360
nr = ceil(r/max(r(:))*360); % normalized rho in degrees
nr = mod(nr + 240 + start,360);  % calibrate angles to stimulus starting position
nr(nr==0) = 1;
if direc 
    nr = 360-nr+1;
end

% image for output
imgR = zeros(2*radius+1, 2*radius+1);
imgG = zeros(2*radius+1, 2*radius+1);
imgB = zeros(2*radius+1, 2*radius+1);

if strcmpi(wheel, 'Polar')
    imgR(r<radius) = cmap(t(r<radius),1);
    imgG(r<radius) = cmap(t(r<radius),2);
    imgB(r<radius) = cmap(t(r<radius),3);
elseif strcmpi(wheel, 'Eccen')
    imgR(r<radius) = cmap(nr(r<radius),1);
    imgG(r<radius) = cmap(nr(r<radius),2);
    imgB(r<radius) = cmap(nr(r<radius),3);
else
    error('What the hell are you asking of me?!');
end

img = imgR;
img(:,:,2) = imgG;
img(:,:,3) = imgB;

% Flip left & right
img = flipdim(img, 2);
