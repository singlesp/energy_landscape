function varargout = ReadPulseSiemensTrio(filename, dcminfo, fs_in, fs_out, varargin)
p = inputParser;
p.addParamValue('pulse_minspace',0.5);
p.addParamValue('trig_offset',0.035);
p.parse(varargin{:});
r = p.Results;

trig_to_peak = r.trig_offset;
pulse_minspace = fix(r.pulse_minspace*fs_out);

%info = dicominfo(dcmfilename);
%info.AcquisitionTime
dvec = sscanf(dcminfo.ContentDate,'%4d%2d%2d');
tvec = sscanf(dcminfo.AcquisitionTime,'%2d%2d%2d.%4d00');
acqtime = datenum([dvec' tvec(1) tvec(2) tvec(3)+tvec(4)/10000]);

fid = fopen(filename,'r');
M = textscan(fid,'%d');
data = M{1}(6:end);

pulse_trig = find(data==5000);
pulse_data = double(data(data < 5000));
pulse_trig = pulse_trig(pulse_trig <= numel(data));
pulse_trig = pulse_trig - 2*(1:numel(pulse_trig))';

%%

pulse_data = interp1(linspace(0,1,numel(pulse_data)),pulse_data,linspace(0,1,numel(pulse_data)*fs_out/fs_in),'spline');
pulse_data = (pulse_data-mean(pulse_data))/std(pulse_data);

pulse_trig = round(pulse_trig * fs_out/fs_in);
pulse_trig = pulse_trig -  round(trig_to_peak*fs_out);

%%
mdhstart = [];
mdhstop = [];
mpcustart = [];
mpcustop = [];

tline = fgetl(fid);
while(ischar(tline))
    words = regexp(tline,' ','split');
    if(strcmpi(words{1},'LogStartMDHTime:'))
        mdhstart = str2num(words{end});
    elseif(strcmpi(words{1},'LogStopMDHTime:'))
        mdhstop = str2num(words{end});
    elseif(strcmpi(words{1},'LogStartMPCUTime:'))
        mpcustart = str2num(words{end});
    elseif(strcmpi(words{1},'LogStopMPCUTime:'))
        mpcustop = str2num(words{end});
    end
    tline = fgetl(fid);
end 
fclose(fid);

pulsdate = 0;
pulstime = 0;
%mdhstart = time of sample 1 (in dicom time units)
logstart = pulsdate + mdhstart/(24*60*60*1000);
logstop = pulsdate + mdhstop/(24*60*60*1000);

%mpcustart = time of sample 1 (in PMU time units)
logstart_pmu = pulsdate + mpcustart/(24*60*60*1000);
logstop_pmu = pulsdate + mpcustop/(24*60*60*1000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%datestr(logstart)
%datestr(logstop)
dur = etime(datevec(logstop),datevec(logstart));
dur_pmu = etime(datevec(logstop_pmu),datevec(logstart_pmu));
%[floor(dur) mod(dur,1)]
%sprintf('%.5f',dur)
%sprintf('%.5f',dur_pmu)

%datestr(acqtime)
%datestr(pulstime)
%datestr(logstart)
%datestr(logstart_pmu)

dv_log = datevec(logstart);
dv_pmu = datevec(logstart_pmu);
dv_acq = datevec(acqtime);
dv_puls = datevec(pulstime);

dv = [dv_log; dv_pmu; dv_acq; dv_puls];

%dv(:,end)-dv(1,end)

tr_from_logstart = dv_acq(end)-dv_log(end);

timevals = [0:numel(pulse_data)-1]/fs_out - tr_from_logstart;
pulse_trig = pulse_trig(pulse_trig >= tr_from_logstart*fs_out);

%%
pulse_sthresh = zeros(size(timevals));
pulse_sthresh(pulse_trig) = 1;

pulse_minspace_filt = filter(ones(pulse_minspace,1),1,pulse_sthresh);

pulse_sthresh((pulse_minspace_filt > 1)) = 0;

pulse_trig = find(pulse_sthresh);

firstpt = find(timevals >= 0);
pulse_trig = pulse_trig - firstpt(1);
physio = pulse_data(timevals >= 0)';

if(nargout == 1)
    varargout = {pulse_trig};
elseif(nargout == 2)
    varargout = {pulse_trig physio};
end