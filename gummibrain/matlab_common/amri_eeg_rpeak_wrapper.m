function varargout = amri_eeg_rpeak_wrapper(EEG,ECG,config,varargin)
% 1. makes sure artifacts are all aligned
% 2. removes Rpeak events in bad epochs
% option a: original ZM on ECG
% option b: pulls out some (usually #1) PC from EEG elec as ECG
% option c: cleans up data by only keeping first few PC, then uses specified
%   EEG electrode as ECG
%
% amri_eeg_rpeak_kj(EEG,'ECG',struct('good_elec',good_elec,'bad_epoch',bad_epoch,'realign',true,'realign_elec_label',realign_label),+display_check);
% amri_eeg_rpeak_kj(EEG,[],struct('method','EEG-PCA','good_elec',good_elec,'bad_epoch',bad_epoch,'realign',true,'realign_elec_label',realign_label),+display_check);
% amri_eeg_rpeak_kj(EEG,[],struct('method','EEG-PCAfilt','pcafilt_num',4,'bcg_electrode',cbc_electrode,'good_elec',good_elec,'bad_epoch',bad_epoch,'realign',true,'realign_elec_label',realign_label),+display_check);

default_config = struct(...
    'method','AMRI', ...
    'good_elec',1:EEG.nbchan, ...
    'bad_epoch',false(1,EEG.pnts), ...
    'bcg_pca',1, ...
    'pcafilt_num', 4, ...
    'bcg_electrode', 'ECG',...
    'realign', true, ...
    'realign_elec_label', 'POz', ...
    'rmarkertype','R', ...
    'canonical_rpeak_to_bcg', 0.21*EEG.srate, ...
    'bcg_filter_size',[], ... 
    'repeat_rpeak',true,...
    'pulse_range',[50 100],...
    'interactive',false);

if(nargin >= 3 && isstruct(config))
    newconfig = default_config;
    fnames = fieldnames(config);
    dfnames = fieldnames(newconfig);
    for f = 1:numel(fnames)
        fcheck = find(strcmpi(fnames{f},dfnames) > 0);
        if(~isempty(fcheck))
            newconfig.(dfnames{fcheck}) = config.(fnames{f});
        end
    end
else
    newconfig = default_config;
end

config = newconfig;

good_elec = config.good_elec;
bad_epoch = config.bad_epoch;
rmarkertype = config.rmarkertype;
bcg_filter_size = config.bcg_filter_size;
canonical_rpeak_to_bcg = config.canonical_rpeak_to_bcg;
use_rmarker = config.repeat_rpeak;
pulse_range = config.pulse_range;

chanlabels = {EEG.chanlocs.labels};

if(isnumeric(config.bcg_electrode))
    config.bcg_electrode = chanlabels{config.bcg_electrode};
end

        
if(strcmpi(config.method,'AMRI'))
    if(~isstruct(ECG))
    
        if(isempty(ECG))
            ECG = config.bcg_electrode;
        end
        
        if(ischar(ECG))
            ECG = {ECG};
        end
        ECG = pop_select(EEG,'channel',ECG);
    end
    %OUTEEG = amri_eeg_rpeak(EEG,ECG,varargin{:});
elseif(strcmpi(config.method,'EEG-PCA'))

    %config.bcg_electrode = 'ECG';
    %bcgidx = StringIndex(chanlabels,config.bcg_electrode);
    bcgidx = [];
    pca_elec = sort(unique([good_elec bcgidx]));
    
    [pc eigvec sv] = runpca(EEG.data(pca_elec,~bad_epoch));
    
    ECG = pop_select(EEG,'channel',pca_elec);
    ECG.data = (eigvec)\EEG.data(pca_elec,:);
    
    for i = 1:numel(ECG.chanlocs)
        ECG.chanlocs(i).labels = num2str(i);
    end
    
    [cbc_elec cbc_elec_rank] = compare_cbc_channels(ECG,1:4,bad_epoch,true,pulse_range);
    
    if(config.bcg_pca <= 0)
        config.bcg_pca = cbc_elec;
    end
    ECG = pop_select(ECG,'channel',config.bcg_pca);


    
elseif(strcmpi(config.method,'EEG-PCAfilt'))

    
    [pc eigvec sv] = runpca(EEG.data(good_elec,~bad_epoch));
    pc_select = eye(size(sv));
    pc_select(:,config.pcafilt_num+1:end) = 0;
    eegfilt = EEG;
    eegfilt.data(good_elec,:) = eigvec*pc_select*((eigvec)\EEG.data(good_elec,:));
    
    %cbc_check_elec = good_elec;
    %cbc_check_elec = good_elec(~cellfun(@isempty,regexpi(chanlabels(good_elec),'^[OP]')));
    %cbc_check_elec = [StringIndex(chanlabels,'ECG') cbc_check_elec];
    %[cbc_elec cbc_elec_rank] = compare_cbc_channels(eegfilt, cbc_check_elec, bad_epoch, true);
    %title(cleantext(sprintf('%s_%04d',subj,scan)));
    %savefig(scan,'pulse_elec_compare',gcf);
    
    ECG = pop_select(eegfilt,'channel',{config.bcg_electrode});
end

ECG = pop_selectevent(ECG,'omittype',rmarkertype,'deleteevents','on');
EEG = pop_selectevent(EEG,'omittype',rmarkertype,'deleteevents','on');

ecgdata_orig = ECG.data;

if(use_rmarker)
    %run once for initial guess at peaks (which just uses a single
    %heartbeat as the template for correlation).  Then average all of the
    %detected beats together to get a cleaner template for peakfinding.
    %Either use the mean heartbeat as template (it will be smooth though
    %and less temporally precise), or find the single initially found peak
    %with the highest correlation to the mean heartbeat and use that as the
    %template.
    
    
    ECG.data(:,bad_epoch) = 0;
    ECG = amri_eeg_rpeak_kj(ECG,ECG,varargin{:},'checking',1,'pulserate',pulse_range);  
    
    %%% remove heartbeats detected during bad data segments
    ev_types = cellfun(@(x)(x(x~=' ')),{ECG.event.type},'UniformOutput',false);
    ev_pulse = StringIndex(ev_types,rmarkertype);
    latency_pulse = [ECG.event(StringIndex(ev_types,rmarkertype)).latency];
    latency_pulse = latency_pulse(latency_pulse > 0 & latency_pulse <= numel(bad_epoch));
    evbad_pulse = ev_pulse(unique(find(bad_epoch(latency_pulse) == 1)));

    
    ECG = pop_selectevent(ECG,'omitevent',evbad_pulse,'deleteevents','on');

    fs = EEG.srate;
    pulse_pre_align = event_related_windows(ECG.data.',eeg_eventfield(ECG,rmarkertype),2*fs,'center');

    bcg_pre_align = event_related_windows(EEG.data.',eeg_eventfield(ECG,rmarkertype),2*fs,'center');
    ECG = realign_cbc_artifact(ECG,1,rmarkertype);
    pulse_post_align = event_related_windows(ECG.data.',eeg_eventfield(ECG,rmarkertype),2*fs,'center');
    
    bcg_post_align = event_related_windows(EEG.data.',eeg_eventfield(ECG,rmarkertype),2*fs,'center');
        
    
    pulse_tval = [0:size(pulse_pre_align,1)-1]/fs - size(pulse_pre_align,1)/2/fs;
    figure;
    ax1=subplotsub(2,3,1,1);
    imagesc(pulse_pre_align);
    ax2=subplotsub(2,3,2,1);
    imagesc(pulse_post_align);
    
    ax3=subplotsub(2,3,1,2);
    plotalpha(pulse_tval,pulse_pre_align,'alpha',.15);
    %hold on;
    %plot(pulse_tval,mean(pulse_pre_align,2),'k','linewidth',3);
    grid on;
    
    ax4=subplotsub(2,3,2,2);
    plotalpha(pulse_tval,pulse_post_align,'alpha',.15);
    %hold on;
    %plot(pulse_tval,mean(pulse_post_align,2),'k','linewidth',3);
    grid on;
    
    ax5=subplotsub(2,3,1,3);
    plot(pulse_tval,mean(pulse_pre_align,2),'r');
    hold on;
    plot(pulse_tval,mean(pulse_post_align,2),'b');
    grid on;
    legend({'pre align','post align'},'location','southoutside');
    
    UniformAxes([ax1 ax2]);
    UniformAxes([ax3 ax4]);
    set([ax3 ax4 ax5],'xlim',[-.5 .5]);
    suptitle('wrapper: rpeaks before/after r-peak realignment');
    OUTEEG = amri_eeg_rpeak_kj(EEG,ECG,'use_rmarker',use_rmarker, varargin{:},'pulserate',pulse_range);  
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% find Rpeaks based on whatever timecourse you provided
    ECG.data(:,bad_epoch) = 0;
    OUTEEG = amri_eeg_rpeak_kj(EEG,ECG,varargin{:},'pulserate',pulse_range);    
    %clear EEG ECG;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% make sure all heartbearts line up against the template (some may have
%%% picked something other than R-peak)
if(config.realign)
    if(isempty(config.realign_elec_label))
        config.realign_elec_label = config.bcg_electrode;
    end
    OUTEEG = realign_cbc_artifact(OUTEEG,config.realign_elec_label,rmarkertype);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% remove heartbeats detected during bad data segments
ev_types = cellfun(@(x)(x(x~=' ')),{OUTEEG.event.type},'UniformOutput',false);
ev_pulse = StringIndex(ev_types,rmarkertype);
latency_pulse = [OUTEEG.event(StringIndex(ev_types,rmarkertype)).latency];
latency_pulse = latency_pulse(latency_pulse > 0 & latency_pulse <= numel(bad_epoch));
evbad_pulse = ev_pulse(unique(find(bad_epoch(latency_pulse) == 1)));


OUTEEG = pop_selectevent(OUTEEG,'omitevent',evbad_pulse,'deleteevents','on');

% %%%%%%%%%%%%%%%%%%%%%%%%%%
ev_types = cellfun(@(x)(x(x~=' ')),{OUTEEG.event.type},'UniformOutput',false);
ev_pulse = StringIndex(ev_types,rmarkertype);
latency_pulse = [OUTEEG.event(ev_pulse).latency];
iri = diff(latency_pulse);
iri_med = round(median(iri));
bcg_window_size = 2*iri_med;
if(isempty(bcg_filter_size) || bcg_filter_size <= 0)
    bcg_filter_size = round(.25*iri_med);
end

%bcg_elec = 1:size(ECG.data,1);
%pulse_windows = event_related_windows(ECG.data',latency_pulse-bcg_window_size/2,bcg_window_size);


bcg_elec = good_elec;
pulse_windows = event_related_windows(OUTEEG.data',latency_pulse-bcg_window_size/2,bcg_window_size);
pulse_window = squeeze(mean(pulse_windows,3));
pulse_window = pulse_window - repmat(mean(pulse_window,1),size(pulse_window,1),1);



bcg_offset_filt = gausswindow(bcg_filter_size);


pulse_gfp = sqrt(sum(pulse_window(:,bcg_elec).^2,2));
pulse_gfp_filt = filtfilt(bcg_offset_filt,1,double(pulse_gfp));

pulse_env = abs(hilbert(pulse_window(:,bcg_elec)));
pulse_env = sqrt(sum(pulse_env.^2,2));
pulse_env_filt = filtfilt(bcg_offset_filt,1,double(pulse_env));

env_window = false(bcg_window_size,1);
env_window(round((-.5*iri_med:.5*iri_med)+bcg_window_size/2)) = true;


%[v maxlat] = max(abs(pulse_gfp_filt(env_window)));
[v maxlat] = max(abs(pulse_env_filt(env_window)));
%[v maxlat] = max(abs(pulse_env(env_window)));
pulse_latency_to_cb_peak = maxlat + min(find(env_window)) - bcg_window_size/2;

if(~(strcmpi(config.method,'AMRI') && strcmpi(config.bcg_electrode,'ECG')))
    pulse_adjust_amount = round(pulse_latency_to_cb_peak - canonical_rpeak_to_bcg);
    for i = 1:numel(ev_pulse)
        ev = ev_pulse(i);
        OUTEEG.event(ev).latency = OUTEEG.event(ev).latency + pulse_adjust_amount;
    end
end
%%%% remove events that have now been shifted out of the valid timerange
ev_types = cellfun(@(x)(x(x~=' ')),{OUTEEG.event.type},'UniformOutput',false);
ev_pulse = StringIndex(ev_types,rmarkertype);
evbad_pulse = ev_pulse(unique(find(latency_pulse < 1 | latency_pulse > size(EEG.data,2))));
OUTEEG = pop_selectevent(OUTEEG,'omitevent',evbad_pulse,'deleteevents','on');

OUTEEG = eeg_checkset(OUTEEG);


bcg_post_align2 = event_related_windows(OUTEEG.data.',eeg_eventfield(OUTEEG,rmarkertype),2*fs,'center');
p_win = size(bcg_pre_align,1)*[1/4 3/4];

%gfp_pow = 2;

gfp_pow = 1;
good_elec = StringIndex({OUTEEG.chanlocs.labels},config.realign_elec_label);

mgfp_pre = squeeze(mean(bcg_pre_align(p_win(1):p_win(2),good_elec,:).^gfp_pow,2));
mgfp_post = squeeze(mean(bcg_post_align(p_win(1):p_win(2),good_elec,:).^gfp_pow,2));
mgfp_post2 = squeeze(mean(bcg_post_align2(p_win(1):p_win(2),good_elec,:).^gfp_pow,2));

[~,mgfp_pre_peak] = max(mgfp_pre,[],1);
[~,mgfp_post_peak] = max(mgfp_post,[],1);
[~,mgfp_post2_peak] = max(mgfp_post2,[],1);

%%
mn = min([mgfp_pre_peak(:); mgfp_post_peak(:); mgfp_post2_peak(:)]);
mx = max([mgfp_pre_peak(:); mgfp_post_peak(:); mgfp_post2_peak(:)]);

xi = linspace(mn,mx,100);
kwidth = 5;
kpre = ksdensity(mgfp_pre_peak,xi,'width',kwidth);
c_pre = diff(sign(diff(kpre))) == -2;
kpost = ksdensity(mgfp_post_peak,xi,'width',kwidth);
c_post = diff(sign(diff(kpost))) == -2;
kpost2 = ksdensity(mgfp_post2_peak,xi,'width',kwidth);
c_post2 = diff(sign(diff(kpost2))) == -2;
% 
% c_pre = (kpre(c_pre) >= 0.25*max(kpre));
% c_post = (kpost(c_post) >= 0.25*max(kpost));
%c_post2 = (kpost2(c_post2) >= 0.25*max(kpost2));

c_pre = xi(c_pre);
c_post = xi(c_post);
c_post2 = xi(c_post2);

d = repmat(mgfp_pre_peak,numel(c_pre),1) - repmat(c_pre(:),1,numel(mgfp_pre_peak));
d_pre = mgfp_pre_peak - c_pre(minindex(abs(d),[],1));
pre_peak_spread = std(d_pre);

d = repmat(mgfp_post_peak,numel(c_post),1) - repmat(c_post(:),1,numel(mgfp_post_peak));
d_post = mgfp_post_peak - c_post(minindex(abs(d),[],1));
post_peak_spread = std(d_post);

d = repmat(mgfp_post2_peak,numel(c_post2),1) - repmat(c_post2(:),1,numel(mgfp_post2_peak));
d_post2 = mgfp_post2_peak - c_post2(minindex(abs(d),[],1));
post2_peak_spread = std(d_post2);
% 
% pre_peak_spread = min(distance(mgfp_pre_peak,c_pre),[],2);
% 
% pre_peak_spread = 1.4826*mad(mgfp_pre_peak*1000/fs,1);
% post_peak_spread = 1.4826*mad(mgfp_post_peak*1000/fs,1);
% post2_peak_spread = 1.4826*mad(mgfp_post2_peak*1000/fs,1);

figure;
subplotsub(3,2,1,1);
plot(xi,kpre*numel(mgfp_pre_peak),'r');
hold on;
hist(mgfp_pre_peak,50);
title(sprintf('orig rpeak.  std=%.4f ms', std(mgfp_pre_peak*1000/fs)));
subplotsub(3,2,2,1);
plot(xi,kpost*numel(mgfp_post_peak),'r');
hold on;
hist(mgfp_post_peak,50);
title(sprintf('align1.  std=%.4f ms',  std(mgfp_post_peak*1000/fs)));
subplotsub(3,2,3,1);
plot(xi,kpost2*numel(mgfp_post2_peak),'r');
hold on;
hist(mgfp_post2_peak,50);
title(sprintf('align2.  std=%.4f ms',  std(mgfp_post2_peak*1000/fs)));
subplotsub(3,2,1,2);
hist(d_pre,50);
title(sprintf('orig rpeak.  std=%.4f ms', std(d_pre*1000/fs)));
subplotsub(3,2,2,2);
hist(d_post,50);
title(sprintf('align1.  std=%.4f ms', std(d_post*1000/fs)));
subplotsub(3,2,3,2);
hist(d_post2,50);
title(sprintf('align2.  std=%.4f ms', std(d_post2*1000/fs)));
suptitle('pulse to MGFP peak');

%%
if(config.interactive)
    latency_pulse = eeg_eventfield(OUTEEG,'R');
    interactive_buffer = 1*OUTEEG.srate;
    
    ecg_window = event_related_windows(ecgdata_orig',latency_pulse,interactive_buffer,'center');
    figure;
    subplot(3,1,1:2);
    himg = imagesc(ecg_window);

    
    subplot(3,1,3);
    himg_zoom = imagesc(ecg_window);
    
    set(himg,'buttondownfcn',{@interactive_click, ecg_window, himg_zoom});
    set(gcf,'windowbuttonmotionfcn',{@(h,e)(interactive_click(himg,e,ecg_window,himg_zoom))});
    
    iri = diff(latency_pulse);
    iri_median = median(iri);
    iri_outliers = iri > 2*iri_median | iri < .5*iri_median;
    iri_outliers_idx = find(iri_outliers);
    fs = ECG.srate;
    figure;
    for i = 1:numel(iri_outliers_idx)
        subplotgrid(numel(iri_outliers_idx),i);
        idx1 = latency_pulse(iri_outliers_idx(i));
        idx2 = latency_pulse(iri_outliers_idx(i)+1);
        
        idxwindow = idx1-interactive_buffer:idx2+interactive_buffer;
        twindow = idxwindow/fs;
        plot(twindow,ecgdata_orig(:,idxwindow)');
        hold on;
        plot(([idx1 idx2])/fs,ecgdata_orig(:,[idx1 idx2])','ro');
        xlim(twindow([1 end]));
    end
    pause;
end

if(nargout == 1)
    varargout = {OUTEEG};
elseif(nargout == 2)
    varargout = {OUTEEG, ECG};
end

function interactive_click(hObj,event,varargin)
ecg = varargin{1};
hzoom = varargin{2};
cp = get(get(hObj,'Parent'),'CurrentPoint');
px = ceil(cp(1,1:2));

zoomsize = 50;
xzoom = px(1) + [-zoomsize/2 zoomsize/2];
set(get(hzoom,'Parent'),'xlim',xzoom);
%pos = getpixelposition(hObj);
%px = px - pos(1:2)

