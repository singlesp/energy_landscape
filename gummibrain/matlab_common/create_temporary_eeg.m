function EEG = create_temporary_eeg(varargin)
if(ischar(varargin{1}))
    s = struct(varargin{:});
else
    s = varargin{1};
end
good_elec = [];
if(isfield(s,'bad_elec'))
    good_elec = 1:numel(s.chanlocs);
    good_elec(s.bad_elec) = [];
    s = rmfield(s,'bad_elec');
end
if(isfield(s,'good_elec'))
    good_elec = s.good_elec;
    s = rmfield(s,'good_elec');
end

if(~isempty(good_elec))
    data = zeros(numel(s.chanlocs),size(s.data,2));
    data(good_elec,:) = s.data;
    s.data = data;
end
    
EEG = eeg_checkset(mergestruct(eeg_emptyset,s));