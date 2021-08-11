function newfile = convert_gz(gzfile, outputdir, prefix, suffix, delete_gzip)

[fdir, fname, fext] = fileparts(gzfile);
if(~regexpmatch(fext,'\.gz$'))
    %warning('This file does not appear to be a gzipped file: %s', gzfile)
    newfile = gzfile;
    return;
end

[~, fname, fext] = fileparts(fname);
newfile_default = [fdir '/' fname fext];

if(~exist('outputdir','var') || isempty(outputdir))
    outputdir = fdir;
end

if(~exist('prefix','var') || isempty(prefix))
    prefix = '';
end

if(~exist('suffix','var') || isempty(suffix))
    suffix = '';
end

if(~exist('delete_gzip','var') || isempty(delete_gzip))
    delete_gzip = false;
end


outfile = gunzip(gzfile,outputdir);

[fdir, fname, fext] = fileparts(outfile{1});
newfile = [fdir '/' prefix fname suffix fext];

%if no prefix/suffix/outdir were given, we need to delete the original gzip
% file so FSL doesn't complain
if(strcmpi(strrep(newfile_default,'\','/'), strrep(newfile,'\','/')))
    if(delete_gzip)
        delete(gzfile);
    end
else
    movefile(outfile,newfile,'f');
end
