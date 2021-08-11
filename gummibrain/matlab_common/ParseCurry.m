function varargout = ParseCurry(filename)
% value_lists = ParseCurry(filename)
% [value_lists property_lists] = ParseCurry(filename)
% creates a nested structure that can be accessed by the section names from
%   the original curry file.
% 
% Ex:
%   value_lists.LOCATION_LIST contains the actual Nx3 location values (if
%       present in this file)
%   property_lists.LOCATION_LIST.LIST_DESCRIPTION returns the LIST_DESCRIPTION property
%       from the LOCATION_LIST section of the curry file
%
% OR (older format):
%
%  [proplist_names proplist_fieldnames proplist_values list_names
%       list_values] = ParseCurry(filename)

fid = fopen(filename,'r'); %open file in read mode

proplist_names = [];
proplist_fieldnames = [];
proplist_values = [];

list_names = [];
list_values = [];

list_order_name = {};
list_order_type = {};

linenum = 0;
while 1
    
    fline = fgetl(fid); %read in 1 line at a time
    linenum = linenum + 1;
    if(~ischar(fline)) %exit if reach end of file
        break;
    end 

    start_match = regexp(fline,'([A-Z0-9_]+) START','tokens');
    start_list_match = regexp(fline,'([A-Z0-9_]+) START_LIST','tokens');

    %%%%% read in a property list
    if(~isempty(start_match) && isempty(start_list_match))
        temp_name = cell2mat(start_match{:});
        list_order_name{end+1} = temp_name;
        list_order_type{end+1} = 'PROPERTIES';
        
        proptemp = textscan(fid,'%s = %[^\n]');
        temp_fieldnames = proptemp{1};
        temp_values = proptemp{2};
        % convert if numeric
        for i = 1:numel(temp_values)
            if(isempty(temp_values{i}))
                continue;
            end
            
            %%% need to use str2num to convert arrays, but it's stupid
            %%% and evaluates so strings like "surface" load a figure
            %%% ... so do this stupid check on the first char
            firstchar = temp_values{i}(1);
            if(findstr(['-' '0':'9'], firstchar))
                numtemp = str2num(temp_values{i});

                if(~isempty(numtemp))
                    temp_values{i} = numtemp;
                end
            end
        end
        proplist_names{end+1} = temp_name;
        proplist_fieldnames{end+1} = temp_fieldnames(1:end-1);
        proplist_values{end+1} = temp_values;
    end

    %%% read in a values list
    if(~isempty(start_list_match))
        temp_name = cell2mat(start_list_match{:});
        list_order_name{end+1} = temp_name;
        list_order_type{end+1} = 'VALUES';
        
        proplist_idx = find(strcmpi(proplist_names,temp_name));
        
        is_numeric = false;
        num_col = 0;
        
        if(~isempty(proplist_idx))
            listtype_idx = find(strcmpi(proplist_fieldnames{proplist_idx},'LIST_TYPE'));
            if(isempty(listtype_idx))
                listtype_idx = find(strcmpi(proplist_fieldnames{proplist_idx},'ListType'));
            end
            listtype = proplist_values{proplist_idx}{listtype_idx};
            
            if(listtype ~= 5)
                is_numeric = true;
                prop_idx = find(strcmpi(proplist_fieldnames{proplist_idx},'LIST_NR_COLUMNS'));
                if(isempty(prop_idx))
                    prop_idx = find(strcmpi(proplist_fieldnames{proplist_idx},'ListNrColumns'));
                end
                num_col = proplist_values{proplist_idx}{prop_idx};
            end
        else
            fpos = ftell(fid);
            fline = fgetl(fid);
            
            m = regexp(fline,'[^0-9eE\-\.\s]','once');
            
            if(isempty(m))
                is_numeric = true;
                num_col = numel(regexp(fline,'[0-9eE\-\.]+'));
            end
            fseek(fid,fpos,'bof');
        end
        
        if(is_numeric)
            if(num_col > 0)
                tmpval = textscan(fid,repmat('%f ',1,num_col));
                tmpval = cellfun(@(x)(x(~isnan(x)&~isinf(x))),tmpval,'uniformoutput',false);
                temp_values = cell2mat(tmpval);
                %temp_values = cell2mat(textscan(fid,repmat('%f ',1,num_col)));

                list_names{end+1} = temp_name;
                list_values{end+1} = temp_values;
            end
        else
            temp_values = [];
            while 1
                fline = fgetl(fid);
                %exit if reach end of list or end of file
                if(~ischar(fline) || ~isempty(findstr(fline,[temp_name ' END_LIST'])))
                    break;
                end 
                temp_values{end+1} = fline;
            end
            
            list_names{end+1} = temp_name;
            list_values{end+1} = temp_values;
        end
    end
end

fclose(fid);

order = struct('name',list_order_name,'type',list_order_type);

if(nargout > 3)
    varargout = {proplist_names proplist_fieldnames proplist_values ...
        list_names list_values};
else

    props = struct();
    for n = 1:numel(proplist_names)
        for f = 1:numel(proplist_fieldnames{n})
            props.(proplist_names{n}).(proplist_fieldnames{n}{f}) = proplist_values{n}{f};
        end
    end

    lists = struct();
    for n = 1:numel(list_names)
        lists.(list_names{n}) = list_values{n};
    end
    
    if(nargout == 1)
        varargout = {lists};
    elseif(nargout == 2)
        varargout = {lists props};
    elseif(nargout == 3)
        varargout = {lists props order};
    end
end