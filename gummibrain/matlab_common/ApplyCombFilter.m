function data = ApplyCombFilter(data, f, fs, fs_filter, db, bw)

if(nargin < 4)
    fs_filter = fs;
end

if(nargin < 5)
    db = -.5;
end

if(nargin < 6)
    bw = .25;
end

comb_f = f;
comb_db = db;
comb_bw = bw;
comb_fs = fs_filter;
comb_o = 8;
CombSpecObj = fdesign.comb('notch','L,BW,GBW,Nsh',round(comb_fs/comb_f),comb_bw,comb_db,comb_o,comb_fs);
CombFilt=design(CombSpecObj);

data = resample(filtfilt(CombFilt.Numerator,CombFilt.Denominator,...
     resample(data,comb_fs,fs)),fs,comb_fs);
 