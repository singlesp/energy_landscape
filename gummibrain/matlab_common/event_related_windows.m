function varargout = event_related_windows(sig, event_idx, windowsize,varargin)
%function [event_windows] = event_related_windows(sig, event_idx, windowsize)
%each column of sig is a timecourse

if(isempty(event_idx))
    if(nargout == 1)
        varargout = {[]};
    elseif(nargout == 2)
        varargout = {[],[]};
    end
    return;
end

if(size(sig,1) == 1)
    sig = sig.';
end

event_idx = event_idx(:);

if(nargin > 3 && (strcmpi(varargin{1},'center') || strcmpi(varargin{1},'centered')))
    event_idx = floor(event_idx - windowsize/2);
end
    
event_idx_invalid = find(event_idx < 1 | event_idx > size(sig,1) - windowsize + 1);
event_idx(event_idx_invalid) = [];

%event_idx = event_idx(event_idx >= 1 & event_idx <= size(sig,1) - windowsize + 1);

if(isempty(event_idx))
    if(nargout == 1)
        varargout = {[]};
    elseif(nargout == 2)
        varargout = {[],[]};
    end
    
    return;
end

window_idx = [0:windowsize-1]';


win_idx = repmat(window_idx,1,numel(event_idx)) + ...
             repmat(event_idx',windowsize,1);
       
event_windows = sig(floor(win_idx(:)),:);
 
event_windows = squeeze(permute(reshape(event_windows,[windowsize size(win_idx,2) size(sig,2) ]),[1 3 2]));
%event_windows = squeeze(reshape(event_windows,[windowsize  size(sig,2) size(win_idx,2)]));

if(nargout == 1)
    varargout = {event_windows};
elseif(nargout == 2)
    varargout = {event_windows, event_idx_invalid};
end