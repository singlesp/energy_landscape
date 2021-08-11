function vals = eeg_eventfield(EEG,event_type,event_field)
if(nargin < 3)
    event_field = 'latency';
end

if(isstruct(EEG) && isfield(EEG,'event'))
    ev = EEG.event;
elseif(isstruct(EEG) && isfield(EEG,'latency'))
    ev = EEG;
end

evtype = cellfun(@(x)(x(x~=' ')),{ev.type},'UniformOutput',false);
evtype_isnum = cellfun(@isnumeric,evtype);
evtype(evtype_isnum) = cellfun(@num2str,evtype(evtype_isnum),'uniformoutput',false);

if(isnumeric(event_type))
    event_type = num2str(event_type);
end

example_val = ev(1).(event_field);
evidx = sort(StringIndex(lower(evtype),lower(event_type)));
if(isnumeric(example_val) || islogical(example_val))
    vals = [ev(evidx).(event_field)];
else
    vals = {ev(evidx).(event_field)};
end
