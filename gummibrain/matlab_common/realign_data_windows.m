function varargout = realign_data_windows(sig, event_latency, windowsize, template_type, varargin)

event_template = [];
if(nargin < 4 || isempty(template_type))
    template_type = 1;
end

if(isnumeric(template_type) && numel(template_type) == windowsize)
    event_template = template_type;
    template_type = [];
end

[event_windows invalid_event] = event_related_windows(sig, event_latency, windowsize,varargin{:});


%event_windows = whitencols(event_windows);
%event_windows(isnan(event_windows)) = 0; %if any pulses were zero'd out (ie: bad epochs)

%event_windows = applyTEO(event_windows,1,EEG.srate,8,40);

valid_latency = event_latency;
valid_latency(invalid_event) = [];

if(isempty(event_template))
    if(isnumeric(template_type) && numel(template_type) == 1)
        event_template = event_windows(:,template_type);        
    elseif(strcmpi(template_type,'mean'))
        event_template = mean(event_windows,2);
    elseif(strcmpi(template_type,'mid'))
        event_template = event_windows(:,round(end/2));
    end
end

xc_all = nan(size(valid_latency));
event_offset = zeros(size(valid_latency));
for i = 1:numel(valid_latency)
    [xc lags] = xcorr(event_windows(:,i),event_template,'coeff');
    if(~any(isnan(xc)))
        [~,ixc] = max(xc);
        event_offset(i) = lags(ixc);
        xc_all(i) = xc(ixc);
    end
    valid_latency(i) = valid_latency(i) + event_offset(i);
end

if(nargout == 1)
    varargout = {valid_latency};
elseif(nargout == 2)
    varargout = {valid_latency,xc_all};
end