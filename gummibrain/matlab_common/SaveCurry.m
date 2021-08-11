function SaveCurry(filename, value_lists, property_lists, list_order)

fid = fopen(filename,'w');
%fid = 1;

for l = 1:numel(list_order)
    listname = list_order(l).name;
    listtype = list_order(l).type;
    
    switch listtype
        case 'PROPERTIES'
            proplist = property_lists.(listname);
            fields = fieldnames(proplist);
            fprintf(fid,'%s START\n',listname);
            for f = 1:numel(fields)
                val = proplist.(fields{f});
                %if(isnumeric(val) && floor(val) == val)
                %    fmt = '%d';
                %elseif(isnumeric(val))
                %    fmt = '%d';
                if(isnumeric(val))
                    fmt = '%d';
                else
                    fmt = '%s';
                end
                fprintf(fid,['   %-20s = ' fmt '\n'],fields{f},val);
            end
            fprintf(fid,'%s END\n\n',listname);
        case 'VALUES'
            fprintf(fid,'%s START_LIST\n',listname);
            vals = value_lists.(listname);
            %if(isnumeric(vals) && all(round(vals(:)) == vals(:)))
            %    fmt = repmat('%d ',1,size(vals,2));
            %    fmt = fmt(1:end-1);
            %elseif(isnumeric(vals))
            %    fmt = repmat('%f ',1,size(vals,2));
            %    fmt = fmt(1:end-1);
            if(isnumeric(vals))
               fmt = repmat('%d ',1,size(vals,2));
               fmt = fmt(1:end-1);
            else
                fmt = '%s';
            end
            if(iscell(vals))
                fprintf(fid,[fmt '\n'],vals{:});
            else
                fprintf(fid,[fmt '\n'],vals);
            end
            fprintf(fid,'%s END_LIST\n\n',listname);
    end
end

if(fid ~= 1)
    fclose(fid);
end