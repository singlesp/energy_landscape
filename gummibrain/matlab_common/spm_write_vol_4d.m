function V = spm_write_vol_4d(V,Vdata)
%
% V = spm_write_vol_4d(V,Vdata)
%
% Write 4D data to .nii file
% V can be either 1x1 or Nx1, where N is size(Vdata,4)
%
% adapted from spm_write_vol.m by KJ

N = size(Vdata,4);

if(numel(V) == N)
elseif(numel(V) == 1)
    V = repmat(V,N,1);
else
    error('V must be 1x1 or same as size(V,4)');
end
    
%%%%%%%%%%%%%%%%
%%% 4D pinfo rescaling modified from spm_write_vol.m
rescal = true;
use_offset = false;

pinfo = V(1).pinfo;
    
if rescal
    % Set scalefactors and offsets
    %-----------------------------------------------------------------------

    
    dt           = V(1).dt(1);
    s            = find(dt == [2 4 8 256 512 768]);
    if isempty(s),
        pinfo(1:2) = [1;0];
    else
        dmnmx        = [0 -2^15 -2^31 -2^7 0 0 ; 2^8-1 2^15-1 2^31-1 2^7-1 2^16-1 2^32-1];
        dmnmx        = dmnmx(:,s);

        mx = max(double(Y(isfinite(Y(:)))));
        mn = min(double(Y(isfinite(Y(:)))));

        if isempty(mx), mx = 0; end;
        if isempty(mn), mn = 0; end;
        if mx~=mn,
            if use_offset,
                pinfo(1,1) = (mx-mn)/(dmnmx(2)-dmnmx(1));
                pinfo(2,1) = (dmnmx(2)*mn-dmnmx(1)*mx)/(dmnmx(2)-dmnmx(1));
            else
                if dmnmx(1)<0,
                    pinfo(1) = max(mx/dmnmx(2),mn/dmnmx(1));
                else
                    pinfo(1) = mx/dmnmx(2);
                end
                pinfo(2) = 0;
            end
        else
            pinfo(1,1) = mx/dmnmx(2);
            pinfo(2,1) = 0;
        end
    end
end

        
%%%%%%%%%%%%%%%%
for i = 1:N
    V(i).n(1) = i;
    V(i).pinfo = pinfo(1:2);
    V(i) = spm_write_vol(V(i),Vdata(:,:,:,i));
end
