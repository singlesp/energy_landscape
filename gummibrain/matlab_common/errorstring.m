function str = errorstring(err,color)
str = {};
str{1} = sprintf('\n%s\n',err.message);
for i = 1:numel(err.stack)
    if(i==1)
        prefix = '>';
    else
        prefix = '';
    end
    st = err.stack(i);
    
    fname = justfilename(st.file,true);
    
    href = sprintf('matlab: opentoline(''%s'',%d,1)',st.file,st.line);
    if(~strcmp(fname,st.name))
        fname = sprintf('%s>%s',fname,st.name);
    end
    txt = sprintf('%s at %d',fname,st.line);
    str{end+1} = sprintf('%-2sIn <a href="%s">%s</a>\n',prefix,href,txt);
end
str = [str{:}];
