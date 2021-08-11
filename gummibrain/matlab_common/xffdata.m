function M = xffdata(varargin)
V = xff(varargin{:});
M = [];
try
    M = V.VoxelData;
catch
end

if(isempty(M))
    try
        M = V.VMRData;
    catch
    end
end

V.ClearObject;