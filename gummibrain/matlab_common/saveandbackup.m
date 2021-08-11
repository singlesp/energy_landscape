function saveandbackup(varargin)
filename=varargin{1};
[d f ext] = fileparts(filename);
filename_backup = fullfile(d,[f '_' timestamp ext]);

tmpstr = strcat('''',varargin(2:end),''', ');
tmpstr = strcat(tmpstr{:});
tmpstr = tmpstr(1:end-1);

str1 = ['save(''' filename ''',' tmpstr ');'];
str2 = ['save(''' filename_backup ''',' tmpstr ');'];

evalin('caller',str1);
evalin('caller',str2);