function [args args_dequote] = parseargs(argstr)
%
% function args = parseargs(argstr)
% function [args dequoted] = parseargs(argstr)
%
% parse argument string and return individual arguments.  arguments
% enclosed in single or double quotes will be preserved as single arguments
%
% optional second output removes quotes from quoted arguments

args = {};
args_dequote = {};

[~,~,~,matchstring,~,~,splitstring] = regexp(argstr,'[\"\''][^\"\'']*[\"\'']');

midx = 1;
for i = 1:numel(splitstring)
    %trim arguments
    s = regexp(splitstring{i},'\s+','split');
    s = s(~cellfun(@isempty,s));
    %s = parsetrim(splitstring{i});
    if(~isempty(s))
        args = [args s];
        args_dequote = [args_dequote s];
    end
    if(midx <= numel(matchstring))
        args = [args matchstring{midx}];
        
        m = regexprep(regexprep(matchstring{midx},'^[\"\''\s]+',''),'[\"\''\s]+$','');
        args_dequote = [args_dequote m];
    end
    midx = midx + 1;    
end
