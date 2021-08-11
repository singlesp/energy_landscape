%% use same filter/downsampling as amri_eeg_gac, for use with data where MRI wasn't running


function [outeeg] = amri_eeg_gac_preproc(eeg,varargin)


if nargin<1
    eval('help amri_eeg_gac_preproc');
    return
end 

%% ************************************************************************
% Defaults
% *************************************************************************

gac_method              =       {'PCA'};    % algorithm name

% trigger markers received from scanner

scanner_trigger.type    =       'slice';    % scanner trigger type (per slice or volume)
scanner_trigger.name    =       {'R128'};   % scanner trigger name 
scanner_trigger.count   =       0;          % total number of scanner triggers
scanner_trigger.index   =       [];         % array of indices to eeg.event
scanner_trigger.latency =       [];         % array of latencies in time points

% scanner events                            
                                            % specify how the artifact would be corrected 
scanner_event.type      =       'slice';    % (by volume or slice)? 
scanner_event.latency   =       [];         % latencies (in time points)
scanner_event.count     =       0;          % total number of scanner events
scanner_event.ignore    =       [];         % list of scanner events to be ignored
scanner_event.tshift    =       [];         % time shift of event latency
scanner_event.epochonset=       -5;         % epoch onset w.r.t. event latency

% fmri acquisition information

fmri.nslice             =       NaN;        % number of slices
fmri.nvolume            =       NaN;        % number of volumes
fmri.trsec              =       NaN;        % volume tr in seconds
fmri.trpnt              =       NaN;        % volume tr in sampling points

% moving window

moving_window.size      =       101;        % size (in number of scanner events)
moving_window.start     =       [];         % index to starting event of each window
moving_window.end       =       [];         % index to ending event of each window
moving_window.index     =       [];         % indices to all events within each window
moving_window.nowin     =       0.3;        % sec, the half window length to exclude neighboring epochs
moving_window.minsize   =       10;         % min size(in number of scanner events)

% parameters for scanner event latency realignment

realign.flag            =       0;          % realign (1) or not (0)
realign.maxtshift       =       1;          % range of time-shift in time points
realign.refscan         =       10;         % scanner event to be realigned to
realign.refchannel      =       1;          % channels used for realignment


% signal processing parameters

upsample_rate           =       1;          % upsample rate for pre-processing
downsample_rate         =       round(eeg.srate/250);          % downsample rate for post-processing

detrend_option          =       NaN;        % no detrend

% taylor series order 

ts.order                =       25;         % order of taylor series expansion

% filters

lowpass.method          =       'ifft';     % filter design method ['ifft'|'firls']
lowpass.filter          =       NaN;        % lowpass filter coefficient
lowpass.cutoff          =       125;        % cutoff frequency (Hz)
lowpass.trans           =       0.04;       % relative transition
lowpass.order           =       NaN;        % filter order

highpass.method         =       'ifft';     % filter design method ['ifft'|'firls']
highpass.filter         =       NaN;        % highpass filter coefficient
highpass.cutoff         =       NaN;        % cutoff frequency (Hz)
highpass.trans          =       0.15;       % relative transition
highpass.order          =       NaN;        % filter order

bandstop.method         =       'ifft';     % filter design method ['ifft'|'butter'|'firls']
bandstop.order          =       NaN;        % order of firls filter
bandstop.trans          =       0.15;       % relative transition (<=0.5)
bandstop.filter         =       NaN;        % filter coefficients
bandstop.width          =       1;          % bandwidth (Hz)
bandstop.lcutoff        =       NaN;        % lowcutoff freq
bandstop.hcutoff        =       NaN;        % highcutoff freq
bandstop.tstart         =       NaN;        % starting time of the time series to be filtered
bandstop.tend           =       NaN;        % ending time of the time series to be filtered

% setting for pca component selection

pca.selection           =       'auto';     % 'auto' | 'fixed'
pca.alpha               =       0.05;       % significance level

% checking intermediate results

checking                =       0;          % checking data
checkfig.position       =       NaN;        % figure position

% flags
flag_verbose            =       0;         % whether print out information

%% ************************************************************************
% Collect keyword-value pairs
% *************************************************************************

for i = 1:2:size(varargin,2) 

    Keyword = varargin{i};
    Value   = varargin{i+1};

    if ~ischar(Keyword)
        printf('amri_eeg_gac(): keywords must be strings'); 
        continue;
    end
    
    if strcmpi(Keyword,'method')
        
        if ischar(Value)
            gac_method{1}=Value;
        elseif iscell(Value)
            gac_method = Value;
        else
            printf('amri_eeg_gac(): invalid value for ''method''');
        end
        
    elseif strcmpi(Keyword,'correctby')
        
        if strcmpi(Value,'slice')
            scanner_event.type='slice';
        else
            scanner_event.type='volume';
        end

    elseif strcmpi(Keyword,'winsize') || strcmpi(Keyword,'win_size')
    
        if isnumeric(Value)
            moving_window.size=round(Value(1));
        elseif isempty(Value)
            moving_window.size=NaN;
        elseif isnan(Value)
            moving_window.size=NaN;
        else
            printf('amri_eeg_gac(): invalid ''winsize'' value');
        end

    elseif strcmpi(Keyword,'trigger.name') || ...
           strcmpi(Keyword,'trigger_name') || ...
           strcmpi(Keyword,'triggername')

        if ischar(Value)
            scanner_trigger.name{1}=Value;
        elseif iscell(Value)
            scanner_trigger.name=Value;
        else
            printf('amri_eeg_gac(): invalid value for ''trigger.name''');
        end

    elseif strcmpi(Keyword,'trigger.type') || ...
           strcmpi(Keyword,'trigger_type') || ...
           strcmpi(Keyword,'triggertype')
       
        if ~ischar(Value)
            printf('amri_eeg_gac(): ''trigger.type'' must be a string');
            return;
        else
            scanner_trigger.type = Value;
        end
        
    elseif strcmpi(Keyword,'ignore')
        scanner_event.ignore = Value;
        
    elseif strcmpi(Keyword,'fmri.nslice') || ...
           strcmpi(Keyword,'nslice') || ...
           strcmpi(Keyword,'numslice') || ...
           strcmpi(Keyword,'nrslice') || ...
           strcmpi(Keyword,'nrofslice') || ...
           strcmpi(Keyword,'nrofslices')

       if ~isempty(Value)
           fmri.nslice=round(Value);
       end
       
    elseif strcmpi(Keyword,'fmri.nvolume') || ...
           strcmpi(Keyword,'nvolume') || ...
           strcmpi(Keyword,'nvol') || ...
           strcmpi(Keyword,'numvolume') || ...
           strcmpi(Keyword,'numvol') || ...
           strcmpi(Keyword,'nrvolume') || ...
           strcmpi(Keyword,'nrvol') || ...
           strcmpi(Keyword,'nrofvolume') || ...
           strcmpi(Keyword,'nrofvol')
       
        if ~isempty(Value)
            fmri.nvolume=round(Value);
        end
        
    elseif strcmpi(Keyword,'fmri.tr') || ...
           strcmpi(Keyword,'tr')
        
       if ~isempty(Value)
            fmri.trsec = Value;
       end

    elseif strcmpi(Keyword,'detrend') || ...
           strcmpi(Keyword,'detrend_opt') || ...
           strcmpi(Keyword,'detrend_option')
       
       detrend_option = Value;
       
    elseif strcmpi(Keyword,'upsample') || ...
           strcmpi(Keyword,'upsample')
        
        if Value<=eeg.srate
            upsample_rate=1;
        else
            upsample_rate = round(Value/eeg.srate);
        end
        
    elseif strcmpi(Keyword,'downsample') || ...
           strcmpi(Keyword,'down_sample')
        
        if Value>=eeg.srate;
            downsample_rate=1;
        else
            downsample_rate = eeg.srate/Value;
        end
        
    elseif strcmpi(Keyword,'lowpass.method')
        
        lowpass.method = Value;
        
    elseif strcmpi(Keyword,'lowpass.trans')
        
        lowpass.trans = Value;
        

    elseif strcmpi(Keyword,'highcutoff') || ...   % compatible
           strcmpi(Keyword,'high_cutoff') || ...  % compatible
           strcmpi(Keyword,'lowpass.cutoff')
        
        if isnumeriic(Value)
            lowpass.cutoff=Value(1);
        else
            lowpass.cutoff = NaN;
        end
        
    elseif strcmpi(Keyword,'highpass.method')
        
        highpass.method = Value;
        
    elseif strcmpi(Keyword,'highpass.trans')
        
        highpass.trans  = Value;
        
    elseif strcmpi(Keyword,'lowcutoff') || ...   % compatible
           strcmpi(Keyword,'low_cutoff') || ...  % compatible
           strcmpi(Keyword,'highpass.cutoff')
        
        if isnumeric(Value)
            highpass.cutoff=Value(1);
        else
            highpass.cutoff=NaN;
        end
        
    elseif strcmpi(Keyword,'realign.flag') || ...
           strcmpi(Keyword,'realign') || ...
           strcmpi(Keyword,'realignment')
        
        if isnumeric(Value)
            realign.flag = Value(1);
        else
            printf('amri_eeg_gac(): ''realign.flag'' must a numeric value');
        end
        
    elseif strcmpi(Keyword,'realign.maxtshift') || ...
           strcmpi(Keyword,'maxtshift') || ...
           strcmpi(Keyword,'max_tshift')
       
        if ~isnumeric(Value)
            printf('amri_eeg_gac(): ''realign.maxtshift'' must be a numeric value');
            return;
        else
            realign.maxtshift=Value;
        end
        
    elseif strcmpi(Keyword,'realign.refscan') || ...
           strcmpi(Keyword,'realign.ref_scan') || ...
           strcmpi(Keyword,'align2scan') || ...
           strcmpi(Keyword,'refscan') || ...
           strcmpi(Keyword,'ref_scan') % Kept for compatibility
       
        if ~isinteger(Value)
            printf('amri_eeg_gac(): ''refscan'' must be a integer value');
            return;
        else
            realign.refscan=Value;    
        end
        
    elseif strcmpi(Keyword,'realign.refchannel') || ...
           strcmpi(Keyword,'realign.refchan') || ...
           strcmpi(Keyword,'realign.ref_channel') || ...
           strcmpi(Keyword,'realign.ref_chan') || ...
           strcmpi(Keyword,'align2chan') || ...
           strcmpi(Keyword,'align2channel') || ...
           strcmpi(Keyword,'refchannel') || ...
           strcmpi(Keyword,'refch') || ...
           strcmpi(Keyword,'ref_channel') || ...
           strcmpi(Keyword,'ref_ch')
       
       if ischar(Value)
           if strcmpi(Value,'all')
               realign.refchannel=1:size(eeg.data,1);
           else
               temp_channel=[];
               for ichan=1:size(eeg.data,1)
                   if strcmpi(Value,eeg.chanlocs(ichan).labels)
                       temp_channel=ichan;
                       break;
                   end
               end
               if isempty(temp_channel)
                   printf(['amri_eeg_gac(): ''refchannel'' ''' Value ''' is not found']);
                   printf(['amri_eeg_gac(): the default refchannel # ' ...
                       int2str(realign.refchannel) ' is used']);
               else
                   realign.refchannel=temp_channel;
               end
           end
           
       elseif isnumeric(Value)
           temp_channel=Value;
           temp_channel((temp_channel>size(eeg.data,1)) | (temp_channel<1))=[];
           if isempty(temp_channel)
               printf('amri_eeg_gac(): the specified ''refchannel'' is not found');
               printf(['amri_eeg_gac(): the default refchannel # ' ...
                       int2str(realign.refchannel) ' is used']);
           else
               realign.refchannel=temp_channel;
           end         
       elseif iscell(Value)
           temp_channel=ones(1,length(Value))*nan;
           for iValue=1:length(Value)
               alabel=Value{iValue};
               if strcmpi(alabel,'all')
                   temp_channel=1:size(eeg.data,1);
                   break;
               end
               for ichan=1:size(eeg.data,1)
                   if strcmpi(alabel,eeg.chanlocs(ichan).labels)
                       temp_channel(iValue)=ichan;
                       break;
                   end
               end
           end
           temp_channel(isnan(temp_channel))=[];
           if isempty(temp_channel)
               printf('amri_eeg_gac(): the specified ''refchannel'' is not found');
               printf(['amri_eeg_gac(): the default refchannel # ' ...
                       int2str(realign.refchannel) ' is used']);
           else
               realign.refchannel=temp_channel;
           end
       else
           printf('amri_eeg_gac(): ''refchannel'' should be either a string or integer');
       end

    elseif strcmpi(Keyword,'check') || strcmpi(Keyword,'checking')
        
        checking = Value;
    elseif strcmpi(Keyword,'verbose')
        if Value>0
            flag_verbose=1;
        else
            flag_verbose=0;
        end
    else
        printf(['amri_eeg_gac(): ' Keyword ' is unknown']);
    end   
end



% if flag_verbose>0
%     printf(['amri_eeg_gac(): moving window of ' ...
%         int2str(moving_window.size) ' ' scanner_event.type 's']);
%     printf(['amri_eeg_gac(): note that a moving window may not be used '...
%             'for some algorithms']);
% end

%% ************************************************************************
%              Filter Setting and design
%
% Using firls for filter design requires signal processing toolbox included
% in the matlab. 
% The default inverse fft method does not require this toolbox. 
% *************************************************************************

% ------------------------
%     lowpass filter
% ------------------------

if isnan(lowpass.cutoff)
    lowpass.cutoff = eeg.srate/downsample_rate/2;
end

if strcmpi(lowpass.method,'firls')
    % design lowpass filter using firls
    minfac = 3;
    min_filtorder = 24;
    fs = eeg.srate;
    F = [0 lowpass.cutoff lowpass.cutoff*(1+lowpass.trans) fs/2];
    F = F * 2 / fs;
    A = [1 1 0 0];
    lowpass.order = max([minfac*fix(fs/lowpass.cutoff) min_filtorder]);
    if mod(lowpass.order,2)==1
        lowpass.order=lowpass.order+1;
    end
    lowpass.filter=firls(lowpass.order,F,A);
elseif strcmpi(lowpass.method,'ifft')
    % design lowpass filter using frequency domain filter

    % sampling frequency
    fs=eeg.srate;
    % number of frequency points
    nfft=2^nextpow2(size(eeg.data,2));
    % frequency vector up to Nygt frequency
    fv = fs/2*linspace(0,1,nfft/2);
    % res in frequency domain
    fres = (fv(end)-fv(1))/(nfft/2-1);
    % frequency-domain filter
    lowpass.filter = ones(1,nfft);
    idxh = round(lowpass.cutoff/fres)+1;
    idxhpt = round(lowpass.cutoff*(1+lowpass.trans)/fres)+1;
    
    idxh = max(min(idxh,nfft),1);
    idxhpt = max(min(idxhpt,nfft),1);
    
    lowpass.filter(idxh:idxhpt)=0.5*(1+sin(pi/2+linspace(0,pi,idxhpt-idxh+1)));
    lowpass.filter(idxhpt:nfft/2)=0;
    lowpass.filter(nfft/2+1:nfft-idxh+1)=lowpass.filter(nfft/2:-1:idxh);
end

% ------------------------
%     highpass filter
% ------------------------

if isnan(highpass.cutoff)
    highpass.cutoff = fmri.nslice/fmri.trsec/2;
end
if strcmpi(highpass.method,'firls')
    % design highpass filter using firls
    minfac=3;
    min_filtorder = 24;
    fs = eeg.srate;
    F = [0 highpass.cutoff*(1-highpass.trans) highpass.cutoff fs/2];
    F = F * 2 / fs;
    A = [0 0 1 1];
    highpass.order = max([minfac*fix(fs/highpass.cutoff) min_filtorder]);
    if mod(highpass.order,2)==1
        highpass.order=highpass.order+1;
    end
    highpass.filter = firls(highpass.order,F,A);
elseif strcmpi(highpass.method,'ifft')
    % design highpass filter for ifft
    % sampling frequency
    fs=eeg.srate;
    % number of frequency points
    nfft=2^nextpow2(size(eeg.data,2));
    % frequency vector up to Nygt frequency
    fv = fs/2*linspace(0,1,nfft/2);
    % res in frequency domain
    fres = (fv(end)-fv(1))/(nfft/2-1);
    % frequency-domain filter
    highpass.filter = ones(1,nfft);
    idxl = round(highpass.cutoff/fres)+1;
    idxlmt = round(highpass.cutoff*(1-highpass.trans)/fres)+1;
    
    idxl = max(min(idxl,nfft),1);
    idxlmt = max(min(idxlmt,nfft),1);
    
    highpass.filter(idxlmt:idxl)=0.5*(1+sin(-pi/2+linspace(0,pi,idxl-idxlmt+1)));
    highpass.filter(1:idxlmt)=0;
    highpass.filter(nfft-idxl+1:nfft)=highpass.filter(idxl:-1:1);
end


%% ************************************************************************
% remove DC for each channel
% *************************************************************************
for ichan=1:size(eeg.data,1)
     eeg.data(ichan,:)=eeg.data(ichan,:)-mean(eeg.data(ichan,:));
end 


%% ************************************************************************
%                       Lowpass Filtering
% *************************************************************************

if lowpass.cutoff<eeg.srate/2
    msg='';
%     if flag_verbose>0
%         fprintf(['amri_eeg_gac(): lowpass filtering (cutoff=' ...
%                 num2str(lowpass.cutoff) 'Hz; ' ...
%                 'order=' int2str(lowpass.order) ') ']);
%     end
    for ichan = 1 : length(eeg.chanlocs)
        if strcmpi(lowpass.method,'ifft') 
             nfft = 2^nextpow2(size(eeg.data,2));
             ts = eeg.data(ichan,:);
             ts_start=1;
             ts_end=length(ts);
             X = fft(ts,nfft);
             X = X.*lowpass.filter;
             ts = real(ifft(X,nfft));
             ts = ts(ts_start:ts_end); % truncate data
             eeg.data(ichan,:)=ts;
        elseif strcmpi(lowpass.method,'firls')
            eeg.data(ichan,:)=filtfilt(lowpass.filter,1,eeg.data(ichan,:));
        end
       
%         if flag_verbose>0
%             if ~isempty(msg)
%                 fprintf(repmat('\b',1,length(msg)));
%             end
%         end
        msg = [int2str(ichan) '/' int2str(size(eeg.data,1))];
%         if flag_verbose>0
%             fprintf(msg);
%         end
    end
%     if flag_verbose>0
%         fprintf('\n');
%     end
end


%% ************************************************************************
%                           Downsample
% *************************************************************************

if downsample_rate>1
    if flag_verbose>0
        printf(['amri_eeg_gac(): downsample by ' int2str(downsample_rate)]);
    end
    eeg.data  = eeg.data(:,1:downsample_rate:size(eeg.data,2));
    eeg.srate = eeg.srate/downsample_rate;
    eeg.pnts  = floor(eeg.pnts/downsample_rate);
    if ~isfield(eeg,'xmin'),eeg.xmin=0; end;
    eeg.xmax = eeg.xmin + (eeg.pnts-1)/eeg.srate;
%     anc.refsig = anc.refsig(1:downsample_rate:length(anc.refsig));
    for ievent=1:length(eeg.event)
        eeg.event(ievent).latency=round(eeg.event(ievent).latency/downsample_rate);
        eeg.event(ievent).latency=max([1 eeg.event(ievent).latency]);
    end
    for ievent=1:length(eeg.urevent)
        eeg.urevent(ievent).latency=round(eeg.urevent(ievent).latency/downsample_rate);
        eeg.urevent(ievent).latency=max([1 eeg.urevent(ievent).latency]);
    end
end


%% ************************************************************************
% Linear Detrend
% *************************************************************************
for ichan=1:size(eeg.data,1)
    eeg.data(ichan,:)=detrend(eeg.data(ichan,:),'linear');
end

%%
outeeg = eeg;
return;





