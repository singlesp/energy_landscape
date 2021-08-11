function cmap=isolum_cbs(n,varargin)
%ISOLUM_CBS Isoluminant-Colormap (printer save & colorblind save)
%
%   isoluminant-Colormap that is "save" for people with colorblindness
%   and save for printing
%
%
%	Written by Matthias Geissbuehler - matthias.geissbuehler@a3.epfl.ch
%	June 2011
%
%   Features:
%     1) All colors have nearly the same luminescence (ideal for lifetime
%        images that will be displayed with an additional transparency map
%        to "mask" places where the lifetime is not well defined)
%     2) Color vision deficient persons can only see reduced color: as much
%        as 10% of adult male persons have a red-green defiency (either
%        Deuteranope  or Protanope) -> as a result they can only distinguish
%        between blue and yellow. A colormap which is "save" for color vision
%        deficient persons is hence only based on these colors.
%        However: people with normal vision DO have a larger space of colors
%        available: it would be a pity to discard this freedom. So the goal
%        must be a colormap that is both using as many colors as possible
%        for normal-sighted people as well as a colormap that will "look"
%        blue-yellow to people with colorblindness without transitions that
%        falsify the information by including a non-distinct transitions
%        (as is the case for many colormaps based on the whole spectrum
%        (ex. rainbow or jet).
%        That's what this colormap here tries to achieve!
%     3) In order to be save for publications, the colormap uses colors that
%        are only from the CMYK colorspace (or at least not too far)
%     4) The map is designed in order to be overlaid with an additional
%        black/white channel (ex. for adding an additional phase-contrast
%        image. However it should be noted that for colorblind people:
%        since blue&yellow are complementary colors there always is a
%        transition in this colormap over gray. This part is difficult to
%        be optimized for colorblind-persons.
%
%
%   Usage:
%   cmap = isolum_cbs(n)
%
%   All arguments are optional:
%
%   n           The number of elements (256)
%   
%   Further on, the following options can be applied
%     'gamma'    The gamma of the monitor to be used (1.8)
%     'isolum'   Isoluminescence property. ('default') or 'strict'
%     'minColor' The absolute minimum value can have a different color
%                ('none'), 'white','black','lightgray', 'darkgray'
%                or any RGB value ex: [0 1 0]
%     'maxColor' The absolute maximum value can have a different color
%     'invert'   (0), 1=invert the whole colormap
%
%   Examples:
%     figure; imagesc(peaks(200));
%     colormap(isolum_cbs)
%     colorbar
%
%     figure; imagesc(peaks(200));
%     colormap(isolum_cbs(256,'gamma',1.8,'minColor','black','maxColor',[0 1 0]))
%     colorbar
%
%     figure; imagesc(peaks(200));
%     colormap(isolum_cbs(256,'invert',1,'minColor','white'))
%     colorbar
%

%   Copyright 2011 Matthias Geissbuehler - matthias.geissbuehler@a3.epfl.ch 
%   $Revision: 1.1 $  $Date: 2011/06/14 12:00:00 $
p=inputParser;
p.addParamValue('gamma',1.8, @(x)x>0);
p.addParamValue('isolum','default', @ischar);
p.addParamValue('minColor','none');
p.addParamValue('maxColor','none');
p.addParamValue('invert',0, @(x)x==0 || x==1);

if nargin==1
    p.addRequired('n', @(x)x>0 && mod(x,1)==0);
    p.parse(n);
elseif nargin>1
    p.addRequired('n', @(x)x>0 && mod(x,1)==0);
    p.parse(n, varargin{:});
else
    p.addParamValue('n',256, @(x)x>0 && mod(x,1)==0);
    p.parse();
end
config = p.Results;
n=config.n;

%the ControlPoints and the spacing between them
switch config.isolum
    case 'strict'
        %the ControlPoints in a very isoluminescence case
        cP(:,1) = [90  190 245]./255; k(1)=1;  %cyan at index 1
        cP(:,2) = [157 157 200]./255; k(2)=13; %purple at index 13
        cP(:,3) = [220 150 130]./255; k(3)=24; %purple at index 24
        cP(:,4) = [245 120 80 ]./255; k(4)=43; %redish at index 43
        cP(:,5) = [180 180 0  ]./255; k(5)=64; %yellow at index 64
        
        for i=1:4  % interpolation between control points, while keeping the luminescence constant
            f{i} = linspace(0,1,(k(i+1)-k(i)+1))';  % linear space between these controlpoints
            ind{i} = linspace(k(i),k(i+1),(k(i+1)-k(i)+1))';
            
            cmap(ind{i},1) = ((1-f{i})*cP(1,i)^config.gamma + f{i}*cP(1,i+1)^config.gamma).^(1/config.gamma);
            cmap(ind{i},2) = ((1-f{i})*cP(2,i)^config.gamma + f{i}*cP(2,i+1)^config.gamma).^(1/config.gamma);
            cmap(ind{i},3) = ((1-f{i})*cP(3,i)^config.gamma + f{i}*cP(3,i+1)^config.gamma).^(1/config.gamma);
        end
    case 'default'
        %the ControlPoints in a bit more colorful variant -> slightly less
        %isoluminescence, but gives a more vivid look
        cP(:,1) = [30  60  150]./255; k(1)=1;  %cyan at index 1
        cP(:,2) = [180 90  155]./255; k(3)=17; %purple at index 17
        cP(:,3) = [230 85  65 ]./255; k(4)=32; %redish at index 32
        cP(:,4) = [220 220 0  ]./255; k(5)=64; %yellow at index 64
        for i=1:3
            f{i}   = linspace(0,1,(k(i+1)-k(i)+1))';  % linear space between these controlpoints
            ind{i} = linspace(k(i),k(i+1),(k(i+1)-k(i)+1))';
        end
        cmap = interp1((1:4),cP',linspace(1,4,64)); % for non-iso points, a normal interpolation gives better results
end

% normal linear interpolation to achieve the required number of points for the colormap
cmap = abs(interp1(linspace(0,1,size(cmap,1)),cmap,linspace(0,1,n)));

if config.invert
    cmap = flipud(cmap);
end

if ischar(config.minColor)
    if ~strcmp(config.minColor,'none')
        switch config.minColor
            case 'white'
                cmap(1,:) = [1 1 1];
            case 'black'
                cmap(1,:) = [0 0 0];
            case 'lightgray'
                cmap(1,:) = [0.8 0.8 0.8];
            case 'darkgray'
                cmap(1,:) = [0.2 0.2 0.2];
        end
    end
else
    cmap(1,:) = config.minColor;
end
if ischar(config.maxColor)
    if ~strcmp(config.maxColor,'none')
        switch config.maxColor
            case 'white'
                cmap(end,:) = [1 1 1];
            case 'black'
                cmap(end,:) = [0 0 0];
            case 'lightgray'
                cmap(end,:) = [0.8 0.8 0.8];
            case 'darkgray'
                cmap(end,:) = [0.2 0.2 0.2];
        end
    end
else
    cmap(end,:) = config.maxColor;
end
