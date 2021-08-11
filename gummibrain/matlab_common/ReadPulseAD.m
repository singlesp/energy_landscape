function varargout = ReadPulseAD(filename, fs_in, fs_out, varargin)
%returns pulse latency relative to first TR
%output is in samples @ fs_out
%
% can also return physio data (pulse timecourse) from first TR to last
    p = inputParser;
    p.addParamValue('pulse_minspace',0.25);
    p.addParamValue('pulse_thresh',3);
    p.parse(varargin{:});
    r = p.Results;

    pulse_mintime = r.pulse_minspace;
    pulse_adthresh = r.pulse_thresh;

    A_stop1 = 12;		% Attenuation in the first stopband
    F_stop1 = 0.5;		% Edge of the stopband
    F_pass1 = 2;	% Edge of the passband
    F_pass2 = 30;	% Closing edge of the passband
    F_stop2 = 35;	% Edge of the second stopband
    A_stop2 = 12;		% Attenuation in the second stopband
    A_pass = 1;		% Amount of ripple allowed in the passband
    BandPassSpecObj = ...
        fdesign.bandpass('Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2', ...
        F_stop1, F_pass1, F_pass2, F_stop2, A_stop1, A_pass, ...
        A_stop2, fs_in);
    BandPassFiltAD = design(BandPassSpecObj,'butter');

    fid = fopen(filename,'r');
    M = textscan(fid,'%f %f %f','headerlines',5);
    fclose(fid);
    
    pulseAD = M{2};
    trAD = M{3};
    trAD = trAD > 4;
    latency_trAD = find(diff(trAD) == 1) + 1;

    curval = 0;
    for i = 1:numel(pulseAD)
        if(isnan(pulseAD(i)))
            pulseAD(i) = curval;
        else
            curval = pulseAD(i);
        end
    end
    
    pulseAD = (pulseAD - mean(pulseAD))/std(pulseAD);
    %pulse_ad = filter(BandPassFiltAD, pulse_ad);
    
    pulseAD = resample(pulseAD,fs_out,fs_in);
    %pulse_ad = filter(BandPassFilt,pulse_ad);
    latency_trAD = floor(latency_trAD * (fs_out/fs_in));
    tr_diff = median(diff(latency_trAD));
    
    %pulse_ad = pulse_ad(latency_tr(1):latency_tr(end));
    pulse_minspace = fix(pulse_mintime*fs_out);
    pulse_overthresh = +(diff(pulseAD > pulse_adthresh) == 1);
    pulse_minspace_filt = filter(ones(pulse_minspace,1),1,pulse_overthresh);
    %plot(pulse_minspace_filt);
    %imagesc([pulse_adthresh pulse_minspace_filt]);
    pulse_overthresh((pulse_minspace_filt > 1)) = 0;
    
    %pulse_trigAD = find(diff(pulseAD > 3) == 1) + 1;
    pulse_trigAD = find(pulse_overthresh);
    %latency_tr = latency_tr - latency_tr(1) + 1;

    latency_pulse = pulse_trigAD - latency_trAD(1);
    latency_pulse = latency_pulse(latency_pulse > 0);
    
    physio = pulseAD(latency_trAD(1):latency_trAD(end)+tr_diff-1);

if(nargout == 1)
    varargout = {latency_pulse};
elseif(nargout == 2)
    varargout = {latency_pulse physio};
elseif(nargout == 3)
    varargout = {latency_pulse physio pulseAD};
end