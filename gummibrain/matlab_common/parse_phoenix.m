function val = parse_phoenix(phx,phx_field)

phx_lines = regexp(phx,'\n','split');
str_line = phx_lines{cellfun(@(s)(~isempty(strfind(s,phx_field))),phx_lines)};
str_val = regexp(str_line,'= (.+)$','tokens');
str_val = str_val{1}{1};

if(~isempty(str2num(str_val)))
    val = str2num(str_val);
else
    val = str_val;
end