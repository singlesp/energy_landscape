function nifti_4d_to_5d(niifile4d,niifile5d)

[niidir,~]=fileparts(niifile4d);
tmpfile=['tmp_' justfilename(niifile4d)];

origdir=pwd;
try
    
    if(~isempty(niidir))
        cd(niidir);
    end

    % Create the nifti file object
    nfd = niftifile(justfilename(niifile4d));
    
    % Open the nifti file for reading
    nfdin = fopen(nfd,'read');
    % nfdin now contains the header information
    
    % initialize the data array
    % Matlab is column major so order the data y,x
    
    [nfdin,data] = fread(nfdin, nfdin.nvox);
    
    % This closes the file, but lets us keep using the header
    nfdin = fclose(nfdin);
    
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create a new file handle with a copy of the other one's header
    nfdout = niftifile(tmpfile,nfdin);
    nfdout.descrip = 'modified';
    
    % Make sure we have the correct data type stored
    nfdout.datatype = class(data);

    nfdout.nt = 1;
    nfdout.nu = 3;
    
    nfdout.du=1;
    % Open for writing
    nfdout = fopen(nfdout,'write'); % at this point the header is written
    
    [nfdout] = fwrite(nfdout, data, nfdout.nvox);
    
    % Close the handle
    nfdout = fclose(nfdout);
catch err
    err
end

cd(origdir);
if(exist(fullfile(niidir,tmpfile),'file'))
    movefile(fullfile(niidir,tmpfile),niifile5d,'f');
end

