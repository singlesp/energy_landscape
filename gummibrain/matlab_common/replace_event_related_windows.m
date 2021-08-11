function varargout = replace_event_related_windows(fullsig, event_windows, event_idx, varargin)
%
%function [fullsig] = replace_event_related_windows(fullsig, event_windows, event_idx, [optional])
%
% event_windows must be consistent shape as return value from event_related_windows:
% if fullsig is 1 column, event_windows must be [eventlen x numevents]
% if fullsig is N columns, event_windows must be [eventlen x N x numevents]
%
% optional:
% replace_event_related_windows(...,'center') means event_idx is the CENTER
%   of the event window, as opposed to the start
% replace_event_related_windows(...,'overlap','linear' (default) or 'pchip')
%   if windows overlap, how to blend the overlapping portions. pchip is
%   simliar to a half cosine.


blendtype = 'linear';
%blendtype = 'pchip';

if(isempty(event_idx))
    if(nargout == 1)
        varargout = {[]};
    elseif(nargout == 2)
        varargout = {[],[]};
    end
    return;
end


windowsize = size(event_windows,1);

event_idx = event_idx(:);

if(nargin > 3)
    if(any(strcmpi(varargin,'center') | strcmpi(varargin,'centered')))
        event_idx = floor(event_idx - windowsize/2);
    end
    blendarg = strcmpi(varargin,'blend') | strcmpi(varargin,'blendtype') | strcmpi(varargin,'overlap');
    if(any(blendarg))
        blendtype = varargin{find(blendarg,1)+1};
    end
end

    

if(size(event_windows,2) == 1)
    event_windows = repmat(event_windows(:),1,numel(event_idx));
end

eventsz = size(event_windows);

if(isempty(fullsig))
    fullsig = zeros(max(event_idx)+windowsize-1,eventsz(end));
end

if(numel(fullsig) <= 3)
    fullsig = zeros(fullsig);
end

if(size(fullsig,1) == 1)
    fullsig = fullsig.';
end

%%%%%%%%%%%%%


event_idx_invalid = find(event_idx < 1 | event_idx > size(fullsig,1) - windowsize + 1);
if(ndims(event_windows) == 3 && size(event_windows,2) == numel(event_idx))
    event_windows(:,event_idx_invalid,:) = [];
elseif(size(event_windows,2) == numel(event_idx))
    event_windows(:,event_idx_invalid) = [];
end
event_idx(event_idx_invalid) = [];


%event_idx = event_idx(event_idx >= 1 & event_idx <= size(fullsig,1) - windowsize + 1);


eventsz = size(event_windows);
sigsz = size(fullsig);

event_nonsingle = max(find(eventsz > 1));
sig_nonsingle = max(find(sigsz > 1));
validsz = (event_nonsingle-sig_nonsingle == 1) || (numel(event_idx) == 1);
if(isempty(event_idx) || ~validsz)
    if(nargout == 1)
        varargout = {[]};
    elseif(nargout == 2)
        varargout = {[],[]};
    end
    
    return;
end

window_idx = [0:windowsize-1]';
win_idx = floor(repmat(window_idx,1,numel(event_idx)) + ...
             repmat(event_idx.',windowsize,1));


dperm = 1:numel(eventsz);

if(numel(event_idx) ~= 1 && event_nonsingle > 2)
    dperm([end-1 end]) = dperm([end end-1]);
end

win_overlap = win_idx(end,1:end-1)-win_idx(1,2:end)+1;
if(any(win_overlap > 0))
    newsig = zeros(size(fullsig));
    newsig_mask = false(size(newsig,1),1);
    for i = 1:size(win_idx,2)
        if(sig_nonsingle > 1)
            ev = event_windows(:,:,i);
        else
            ev = event_windows(:,i);
        end
        
        win1 = win_idx(:,i);
        exmask = true(size(win1));
        if(i > 1)
            exmask = exmask & (win1 > win_idx(end,i-1));
        end
        if(i < size(win_idx,2))
            exmask = exmask & (win1 < win_idx(1,i+1));
        end
        
        winmask = ones(size(win1));
        presz = find(exmask,1) - 1;
        postsz = numel(winmask)-find(exmask,1,'last');
        
        if(isempty(presz))
            presz = 0;
        end
        if(isempty(postsz))
            postsz = 0;
        end        
        
        preblend = interp1([-1 0 1 2],[0 0 1 1],linspace(0,1,presz),blendtype);
        postblend = interp1([-1 0 1 2],[1 1 0 0],linspace(0,1,postsz),blendtype);
        
        winmask(1:presz) = preblend;
        winmask(end-postsz+1:end) = postblend;

        %winmask = [linspace(0,1,min(find(exmask))-1) ones(1,sum(exmask)) linspace(1,0,numel(exmask)-max(find(exmask)))];
        ev = bsxfun(@times,ev,winmask(:));
        
        newsig(win1,:) = newsig(win1,:) + reshape(permute(ev,dperm),[],sigsz(end));
        newsig_mask(win1) = true;
    end

    fullsig(newsig_mask,:) = newsig(newsig_mask,:);
else
    fullsig(win_idx(:),:) = reshape(permute(event_windows,dperm),[],sigsz(end));
end
%event_windows = squeeze(permute(reshape(event_windows,[windowsize size(win_idx,2) size(fullsig,2) ]),[1 3 2]));
%event_windows = squeeze(reshape(event_windows,[windowsize  size(fullsig,2) size(win_idx,2)]));

if(nargout == 1)
    varargout = {fullsig};
elseif(nargout == 2)
    varargout = {fullsig, event_idx_invalid};
end