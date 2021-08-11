%A = xml2struct('retinfo.txt');
function S = xml2struct_kj(xml)

if(ischar(xml) && exist(xml,'file'))
    xml = xml2struct(xml);
elseif(isstruct(xml))
else
    disp(xml);
    error('unknown input format');
end
S = [];
fields = fieldnames(xml);


for f = 1:numel(fields)
    name = fields{f};
    val=xml.(fields{f});
    if(isstruct(val))
        vfields = fieldnames(val);
        if(numel(vfields)==1 && strcmp(vfields{1},'Text'))
            val = val.Text;
        else
            val = xml2struct_kj(val);
        end
    else
        
    end
    if(isempty(val))
        val = [];
    elseif(ischar(val))
        num = str2num(val); %#ok<ST2NM>
        if(~isempty(num))
            val = num;
        end
    end
    S.(name) = val;
end
