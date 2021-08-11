function casenames = ListSwitchCases(mfile)


if(isunix)
    cmdprefix = '';
else
    cmdprefix = 'c:\cygwin\bin\';
end

perlstr = 'case\s''';
%sedstr = 's/^\s*case\s//';
cmdstr = sprintf('%scat %s | %sgrep --perl-regexp \"%s\"',cmdprefix,mfile,cmdprefix,perlstr);
[result,output] = system(cmdstr);
if(result == 0)
    S = regexp(output,'[\r\n]+','split');
    S = flatten(regexpi(S,'^\s*case\s\''?([^\'']+)\''?\s?$','tokens'));
    bnum = cellfunc(@isnumeric,S);
    S(bnum) = cellfunc(@str2num,S(bnum));
    casenames = S;
else
    casenames = {};
    fid = fopen(mfile,'r');
    str = fgetl(fid);
    while ischar(str)
        S = flatten(regexpi(str,'\s*case\s\''?([^\'']+)\''?\s?$','tokens'));
        if(~isempty(S))
            num = str2num(S{1});
            if(~isempty(num))
                casenames{end+1} = num;
            else
                casenames{end+1} = S{1};
            end
        end
        str = fgetl(fid);
    end
    fclose(fid);
end
