function data = eeg_interp_subj(data,varargin)
% data can be either EEGlab data structure or a TxChannels array
% eeg_interp_subj(data_array, 'subject', <subj>) 
% eeg_interp_subj(EEGstruct, 'subject', <subj>) 
% eeg_interp_subj(<data>, 'interp_elec', electrodes,'chanlocs',<chanlocs>) 

p = inputParser;
p.addParamValue('subject','');
p.addParamValue('interp_elec',[]);
p.addParamValue('exclude_elec',[]);
p.addParamValue('chanlocs',[]);
p.addParamValue('scan',[]);
p.addParamValue('method','spherical');

p.parse(varargin{:});
r = p.Results;

subj = r.subject;
arg_interp_elec = r.interp_elec;
arg_exclude_elec = r.exclude_elec;
arg_chanlocs = r.chanlocs;
arg_scan = r.scan;
interp_method = r.method;

chanlocs = [];
chanlabels = [];
exclude_elec = [];
%good_elec = [];
interp_elec = [];

if(iscell(subj) && numel(subj) == 1)
    subj = subj{1};
end

if(~isempty(subj))
    expinfo = rivparams(subj);
    matfile = sprintf('%s/%s/%s_%04d_filtered%s.mat',rivsimdir,subj,subj,expinfo.scans(1),expinfo.suffix);
    load(matfile,'chanlocs');
    chanlabels = {chanlocs.labels};
    [exclude_elec, ~, interp_elec] = GetBadElectrodesBA64(subj,chanlabels,true);
end

if(~isempty(arg_scan))
    sc = find(expinfo.scans == arg_scan);
    if(isempty(sc))
        sc = 1;
    end
else
    sc = 1;
end

if(~isempty(arg_interp_elec))
    interp_elec = arg_interp_elec;
end

if(~isempty(arg_exclude_elec))
    exclude_elec = arg_exclude_elec;
end

if(~isempty(arg_chanlocs))
    chanlocs = arg_chanlocs;
    chanlabels = {chanlocs.labels};
end

if(isempty(chanlocs) || (isempty(exclude_elec) && isempty(interp_elec)))
    error('missing inputs');
end

if(iscell(interp_elec) && isempty(interp_elec))
    interp_elec = [];
end

if(iscell(exclude_elec) && isempty(exclude_elec))
    exclude_elec = [];
end

if(iscell(interp_elec) && (iscell(interp_elec{1}) || isnumeric(interp_elec{1})))
    interp_elec = interp_elec{sc};
end

if(iscell(interp_elec) || ischar(interp_elec))
    interp_elec = StringIndex(chanlabels,interp_elec);
end

if(iscell(exclude_elec) || ischar(exclude_elec))
    exclude_elec = StringIndex(chanlabels,exclude_elec);
end

datasize = size(data);
if(~isempty(interp_elec))
    if(isstruct(data))
        EEG = data;
        orig_chanlocs = EEG.chanlocs;
    else
        if(any(size(data) == numel(data)))
            data = data(:);
        else
            data = data.';
        end
        EEG = eeg_emptyset;
        EEG.data = data;
        EEG.nbchan = size(EEG.data,1);
        EEG.pnts = size(EEG.data,2);
        EEG.trials = 1;
    end
    dataclass = class(EEG.data);
    EEG.chanlocs = chanlocs;

    tmp = EEG.data(exclude_elec,:);
    EEG = eeg_interp(EEG,unique([exclude_elec interp_elec]),interp_method);
    EEG.data(exclude_elec,:) = tmp;

    if(strcmp(dataclass,'double'))
        EEG.data = double(EEG.data);
    end

    if(isstruct(data))
        EEG.chanlocs = orig_chanlocs;
        data = EEG;
    else
        data = EEG.data.';
        if(any(size(data) == numel(data)))
            data = reshape(data,datasize);
        end
    end
end