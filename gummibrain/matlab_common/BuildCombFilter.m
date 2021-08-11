function combfilt = BuildCombFilter(f, fs, fs_filter)

comb_db = -.5;
comb_bw = .25;
comb_fs = fs_filter;
comb_o = 8;
CombFilt = {};
filtstr = 'resample(x,comb_fs,fs)';
%L = max(round(comb_fs./f));
%filtstr = 'x';
for i = 1:numel(f)
    comb_f = f(i);
    L = round(comb_fs/comb_f);
    CombSpecObj = fdesign.comb('notch','L,BW,GBW,Nsh',L,comb_bw,comb_db,comb_o,comb_fs);
    CombFilt{i}=design(CombSpecObj);
    filtstr = sprintf('filtfilt(CombFilt{%d}.Numerator,CombFilt{%d}.Denominator,%s)',i,i,filtstr);
end
filtstr = sprintf('@(x)(resample(%s,fs,comb_fs));',filtstr);
%filtstr = sprintf('@(x)(%s);',filtstr);

combfilt = eval(filtstr);