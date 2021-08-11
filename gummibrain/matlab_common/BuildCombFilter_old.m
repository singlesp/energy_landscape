function combfilt = BuildCombFilter(f, fs, fs_filter)

comb_f = f;
comb_db = -.5;
comb_bw = .25;
comb_fs = fs_filter;
comb_o = 8;
CombSpecObj = fdesign.comb('notch','L,BW,GBW,Nsh',round(comb_fs/comb_f),comb_bw,comb_db,comb_o,comb_fs);
CombFilt=design(CombSpecObj);
combfilt = @(x)(resample(filtfilt(CombFilt.Numerator,CombFilt.Denominator,...
    resample(x,comb_fs,fs)),fs,comb_fs));
 