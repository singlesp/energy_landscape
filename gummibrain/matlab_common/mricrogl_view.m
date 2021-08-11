function varargout = mricrogl_view(varargin)
p = inputParser;

p.addRequired('func',@(x)exist(x,'file')>0);
p.addRequired('anat',@(x)exist(x,'file')>0);
p.addParamValue('azel',[60 30]);
p.addParamValue('mip',false);
p.addParamValue('wait',false);
p.addParamValue('display',true);
p.addParamValue('thresh',.0147);
p.addParamValue('scriptfile',[]);

p.parse(varargin{:});
r = p.Results;

anatomical_file = r.anat;
activation_file = r.func;
viewazel = r.azel;
do_wait = r.wait;
pthresh = r.thresh;
do_mip = r.mip;
scriptfile = r.scriptfile;
do_display = r.display;

if(do_mip)
    mipstr = 'true';
else
    mipstr = 'false';
end

M = load_nii(activation_file);

if(pthresh > 0)
    histsize = 1000;
    mvals = M.img(:);
    mvals = mvals(mvals ~= 0);

    [hc x] = hist(mvals,histsize);
    hc_pdf = hc/sum(hc);
    hc_cdf = cumsum(hc_pdf);

    [~,xminidx] = min(abs(hc_cdf-pthresh));
    [~,xmaxidx] = min(abs(1-hc_cdf-pthresh));
    maxval = max(abs(x([xminidx xmaxidx])));
else
    maxval = max(abs(M.img(:)));
end

maxval = maxval*1.1;


template_dict = {
    '%%%%ANATFILE%%%%',anatomical_file
    '%%%%SPMFILE%%%%',activation_file
    '%%%%MAXVAL%%%%', sprintf('%f',maxval)
    '%%%%COLORMAP%%%%','plusminus2'
    '%%%%AZIMUTH%%%%',sprintf('%d',round(viewazel(1)))
    '%%%%ELEVATION%%%%',sprintf('%d',round(viewazel(2)))
    '%%%%OVERLAYFORMVISIBLE%%%%','false'
    '%%%%MAXPROJECTION%%%%',mipstr
    };

mricrogl_exe = cleantext('C:\Program Files (x86)\mricrogl\mricrogl.exe');
glstemplate = 'mricrogl_template.gls';
if(~isempty(scriptfile))
    glsoutput = scriptfile;
else
    %glsoutput = fullfile(tempdir,sprintf('tmp_mricrogl_script_%s.gls',datestr(now,'YYYYmmdd-HHMM')));
    %glsoutput = sprintf('%s_mricrogl_script.gls',tempname);
    [dn,fn,ext] = fileparts(activation_file);
    glsoutput = sprintf('%s/%s_mricrogl_view.gls',dn,fn);
end

glstext = {};
fid = fopen(glstemplate,'r');
tline = fgetl(fid);
while(ischar(tline))
    glstext{end+1} = tline;
    tline = fgetl(fid);
end
fclose(fid);

template_dict = fixslash(template_dict);
for i = 1:size(template_dict,1)
    glstext = regexprep(glstext,template_dict{i,1},template_dict{i,2});
end

fid = fopen(glsoutput,'w');
fprintf(fid,'%s\n',glstext{:});
fclose(fid);

if(do_display && exist(mricrogl_exe,'file'))
    if(do_wait)
        system(['"' mricrogl_exe '" "' glsoutput '"']);
    else
        system(['start "' mricrogl_exe '" "' glsoutput '"']);
    end
end

if(nargout > 0)
    varargout = {glsoutput};
end

