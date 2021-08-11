function [dipole model EEGOUT fieldtrip_extras] = dipfit_eeglab(DATA,chanlocs,gridspacing, varargin)

if(nargin < 3 || isempty(gridspacing))
    gridspacing = 7;
end

ncomps = size(DATA,2);

% dummy EEG dataset
% ---------------------
EEG          = eeg_emptyset;
EEG.data     = rand(size(DATA,1), 1000);
EEG.srate    = 1000;
EEG.nbchan   = size(DATA,1);
EEG.pnts     = 1000;
EEG.trials   = 1;
EEG.chanlocs = chanlocs;
EEG.icawinv    = [ DATA DATA ];
EEG.icaweights = zeros(size([ DATA DATA ]))';
EEG.icasphere  = zeros(size(DATA,1), size(DATA,1));

EEG.icaact     = EEG.icaweights*EEG.icasphere*EEG.data(:,:);
EEG.icaact     = reshape( EEG.icaact, size(EEG.icaact,1), size(EEG.data,2), size(EEG.data,3));
%EEG.dipfit = [];
EEG            = eeg_checkset(EEG);

dipsettings = {...
    'hdmfile','standard_BESA.mat',...
    'coordformat','Spherical',...
    'mrifile','avg152t1.mat',...
    'chanfile','standard-10-5-cap385.elp',...
    'coord_transform',[0 0 0 0 0 0 85 85 85] ,...
    'chansel',[1:8]};

EEG = pop_dipfit_settings( EEG, dipsettings{:});

dipfitdefs;

xg = -floor(meanradius):gridspacing:floor(meanradius);
yg = -floor(meanradius):gridspacing:floor(meanradius);
zg = 0:gridspacing:floor(meanradius);

EEG = pop_dipfit_gridsearch( EEG, [1:ncomps], xg, yg, zg, 100);

dipole = EEG.dipfit.model;
model = EEG.dipfit;

if(nargout > 2)
    EEGOUT = EEG;
end

if(nargout > 3)
    fieldtrip_extras = EEG.fieldtrip;
end