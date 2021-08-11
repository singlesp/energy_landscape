function varargout = unixvmpath(pathstr, dirmaps_new, reverse)

if(exist('reverse','var') && (ischar(reverse) || reverse > 0))
    reverse = true;
else
    reverse = false;
end

if(exist('dirmaps_new','var'))
    dirmaps = dirmaps_new;
else
    dirmaps = unixvmpath_defaults;
end


if(~isempty(pathstr))
    dirmaps_struct = struct('local',dirmaps(1:2:end),'remote',dirmaps(2:2:end));

    pathstr = strrep(pathstr,'\','/');
    
    pathstr = strrep(pathstr,'(','\(');
    pathstr = strrep(pathstr,')','\)');
    %VirtualBox mounts the D:/ drive to HostData, so swap d:/data for
    %directory before calling
    for d = 1:numel(dirmaps_struct)
        
        if(reverse)
            %convert remote paths back to local
            dtmp = dirmaps_struct(d).local;
            if(iscell(dtmp))
                dtmp = dtmp{1};
            end
            pathstr = strrep(pathstr,dirmaps_struct(d).remote,dtmp);            
        else
            %convert local paths to remote
            dtmp = dirmaps_struct(d).local;
            if(ischar(dtmp))
                dtmp = {dtmp};
            end
            for i = 1:numel(dtmp)
                pathstr = strrep(pathstr,dtmp{i},dirmaps_struct(d).remote);
            end
        end
    end
end

if(nargout == 2)
    varargout = {pathstr, dirmaps};
else
    varargout = {pathstr};
end
