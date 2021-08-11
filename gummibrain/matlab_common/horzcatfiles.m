function varargout = horzcatfiles(outputfile,inputfiles,numlines)

if(nargin < 3 || isempty(numlines))
    numlines = [];
end

if(ischar(inputfiles))
    inputfiles = {inputfiles};
end

txt_all = [];
ncols_all = [];

for f = 1:numel(inputfiles)
    
    [inputfile_dir,inputfile_name,inputfile_ext] = fileparts(inputfiles{f});
    fid = fopen(inputfiles{f},'r');

    tline = fgetl(fid);

    txt = [];
    while(ischar(tline))
        txt{end+1} = tline;
        tline = fgetl(fid);
    end
    fclose(fid);
    
    txt = txt(:);
    if(isempty(numlines))
        numlines = numel(txt);
    end
    
    [~,ncols] = sscanf(txt{1},'%s ');
    
    if(numel(txt) > numlines)
        warning(sprintf('file too long. truncating. lines=%d, target=%d :\n\t%s/\n\t\t%s',numel(txt),numlines,inputfile_dir,[inputfile_name inputfile_ext]));
        txt = txt(1:numlines);
    elseif(numel(txt) < numlines)
        warning(sprintf('file too short. zero-padding. lines=%d, target=%d :\n\t%s/\n\t\t%s',numel(txt),numlines,inputfile_dir,[inputfile_name inputfile_ext]));
        zline = sprintf(repmat('%f ',1,ncols),zeros(ncols,1));
        padtxt = repmat({zline},numlines-numel(txt),1);

        txt = [txt(:); padtxt(:)];
    end

    if(isempty(txt_all))
        txt_all = txt;
    else
        txt_all = strcat(cellfun(@(x)([x ' ']),txt_all,'uniformoutput',false),txt);
    end
    ncols_all = [ncols_all ncols];
end

fid = fopen(outputfile,'w');
fprintf(fid,'%s\n',txt_all{:});
fclose(fid);

if(nargout > 0)
    varargout = {ncols_all};
end