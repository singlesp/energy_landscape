function varargout = savebackup(varargin)
[d f ext] = fileparts(varargin{1});
filename = fullfile(d,[f '_' timestamp ext]);

tmpstr = strcat('''',varargin(2:end),''', ');
tmpstr = strcat(tmpstr{:});
tmpstr = tmpstr(1:end-1);

str1 = ['save(''' filename ''',' tmpstr ');'];
evalin('caller',str1);

if(nargout > 0)
    varargout = {filename};
end