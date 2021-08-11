function EEG = eeg_interp_exclude(EEG, interp_elec, bad_elec, varargin)

tmp = EEG.data(bad_elec,:);
EEG = eeg_interp(EEG,unique([bad_elec interp_elec]),varargin{:});
EEG.data(bad_elec,:) = tmp;