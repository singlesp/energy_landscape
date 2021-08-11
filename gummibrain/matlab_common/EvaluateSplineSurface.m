function [P Smat Sinv] = EvaluateSplineSurface(values,varargin)
%P = EvaluateSplineSurface(values,splinefile)
%P = EvaluateSplineSurface(values,G,gx)
%P = EvaluateSplineSurface(values,S)
%[P S] = EvaluateSplineSurface([],splinefile) - Returns spline matrix
%
%values = [num_elec x num_timepoint] electrode values (column = timepoint)
%splinefile = spline .mat file created by headplot('setup')
%
%    or 
%
%S = spline matrix
%
%    or
%
%G is a [num_elec x num_elec] matrix, sort of a spatial blur
%gx is a [num_surf_vert x num_elec] matrix that does actual spline interp
%
%Both of these are contained within the spline .mat file created by
%headplot('setup')
%
%code was extracted/modified from EEGLAB headplot function

S = [];
include_elec = [];
G = [];
gx = [];

if(nargin == 2 && ischar(varargin{1}))
    splinefile = varargin{1};
    load(splinefile); %loads in, among other things, G and gx
elseif(nargin == 2)
    if(isstruct(varargin{1}))
        G = varargin{1}.G;
        gx = varargin{1}.gx;
    else
        S = varargin{1};
    end
elseif(nargin == 3)
    if(isstruct(varargin{1}))
        G = varargin{1}.G;
        gx = varargin{1}.gx;        
        include_elec = varargin{2};
    else
        G = varargin{1};
        gx = varargin{2};
    end
elseif(nargin == 4)
    G = varargin{1};
    gx = varargin{2};
    include_elec = varargin{3};
else
    error('invalid inputs.  see help file for details');
end
    
if(~isempty(values) && isempty(include_elec))
    %exclude channels that are entirely NaN (alternate way to specify included channels)
    include_elec = find(all(~isnan(values),2));
end

if(isempty(include_elec))
    include_elec = 1:size(G,1);
end

if(~isempty(include_elec) && ~isempty(gx))
    fullsz = size(gx);
    G = G(include_elec,include_elec);
    gx = gx(:,include_elec);
end

if(isempty(S)) %loading spline from file or given G,gx

    enum = size(G,1);
    lamd = 0.1;
    S = gx*pinv([(G + lamd);ones(1,enum)]); %spline matrix
    
    S = S(:,1:end-1); %last column was just to help invert
    
    if(~isempty(include_elec))
        Stmp = zeros(fullsz);
        Stmp(:,include_elec) = S;
        S = Stmp;
    end
end

if(isempty(values)) %might only need to return Smat
    P = [];
else
    [sr sc] = size(S);
    [r c] = size(values);
    if(sc == c && sc ~= r)
        values = values.';
    end    
    meanval = nanmean(values,1);
    values = values - ones(size(values,1),1)*meanval; % make mean zero
    %P = S * [values(:);0] + meanval; % fixing division error
    
    P = S(:,include_elec) * values(include_elec,:) + ones(size(S,1),1)*meanval;
end

if(nargout > 1)
    Smat = S;
end

if(nargout > 2)
    Sinv = G*pinv([gx+lamd;ones(1,size(gx,2))]);
    Sinv(:,end) = [];
    if(~isempty(include_elec))
        Stmp = zeros(fullsz(2),fullsz(1));
        Stmp(include_elec,:) = Sinv;
        Sinv = Stmp;
    end
end
