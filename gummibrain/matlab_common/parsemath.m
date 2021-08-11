function [terms friendlyname] = parsemath(str)

str = strrep(str,'((minus))',' - ');
str = strrep(str,'((plus))',' + ');
str = strrep(str,'((times))','*');

segs = regexp(str,' ','split');
segs = segs(~cellfun(@isempty,segs));

is_op = ~cellfun(@isempty,regexp(segs,'^[-+]$'));
is_term = ~is_op;

curop = '+';

expression_str = '';
terms = struct();
termcount = 0;
termname = [];
for s = 1:numel(segs)
    if(is_op(s))
        curop = segs{s};
        continue;
    end
    termcount = termcount+1;
    
    term = segs{s};
    scaleop_lhs = regexp(term,'(^[\-]?[0-9]+)([\*/])([a-zA-Z_\.][a-zA-Z0-9_\.]*$)','tokens');
    scaleop_rhs = regexp(term,'(^[a-zA-Z_\.][a-zA-Z0-9_\.]*)([\*/])([\-]?[0-9]+$)','tokens');
    
    scalar = [];
    varname = [];
    scaleop = [];
    if(~isempty(scaleop_lhs))
        scalar = str2num(scaleop_lhs{1}{1});
        varname = scaleop_lhs{1}{3};
        scaleop = scaleop_lhs{1}{2};
    elseif(~isempty(scaleop_rhs))
        scalar = str2num(scaleop_rhs{1}{3});
        varname = scaleop_rhs{1}{1};
        scaleop = scaleop_rhs{1}{2};
    end
    
    if(isempty(scalar))
        scalar = 1;
        scaleop = '*';
        varname = term;
    end
  
    if(scaleop == '*')
        scalar = scalar;
    elseif(scaleop == '/')
        scalar = 1/scalar;
    end
    
    opstr = '';
    if(curop == '+')
        opstr = '((plus))';
    elseif(curop == '-')
        scalar = -scalar;
        opstr = '((minus))';
    end
    
    if(~(abs(scalar)==1))
        opstr = sprintf('%s%d((times))',opstr,abs(scalar));
    end

    if(termcount == 1 && abs(scalar) == 1 && curop == '+')
        opstr = '';
    end
    terms(termcount).name = varname;
    terms(termcount).scalar = scalar;
    
    expression_str = sprintf('%s%s%s',expression_str,opstr,varname);
end

if(nargout > 1)
    friendlyname = expression_str;
end