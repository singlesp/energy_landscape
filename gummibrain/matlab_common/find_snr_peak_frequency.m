function varargout = find_snr_peak_frequency(Fsnr,freqtags,fs,searchsizehz,freqbinsize,override_thresholding)

if(nargin < 6)
    override_thresholding = false;
end

if(~isreal(Fsnr))
    Fsnr = abs(Fsnr);
end
if(size(Fsnr,1) == 1)
    Fsnr = Fsnr.';
end

if(size(Fsnr,2) > 1)
    fprintf('\nfinding peak frequency in multichannel mode (best freq across all columns)\n');
end

if(nargin < 5 || isempty(freqbinsize))
    freqbinsize = fs/size(Fsnr,1);
end
searchsize = ceil(.5 * searchsizehz / freqbinsize);


fidx_orig = floor(freqtags/freqbinsize);

fidx = fidx_orig;
chidx = -1*ones(size(fidx_orig));

for f = 1:numel(fidx)
	tmpidx = fidx(f) + [-searchsize:searchsize];
    if(numel(tmpidx) < 2)
        fidx(f) = fidx_orig(f);
    end
    if(tmpidx(1) < 1 || tmpidx(end) > size(Fsnr,1))
        continue;
    end
	tmpsnr = Fsnr(tmpidx,:);
	tmpsnr_norm = tmpsnr./abs(max(tmpsnr(:)));

	[tmpsnr_norm sortidx] = sort(tmpsnr_norm,'descend');
	
    [~,chansortidx] = sort(tmpsnr_norm(1,:),'descend');
    
    %try to find column (channel) where there is a VERY large freq with the
    %rest all smaller
    idx = find(all(tmpsnr_norm(2:end,chansortidx) < .75,1),1,'first');
    if(isempty(idx))
        idx = find(all(tmpsnr_norm(3:end,chansortidx) < .75,1),1,'first');
    end
    %if these are both empty, there isn't a strong peak within this window,
    %so return the original frequency (unless overriding)
    if(isempty(idx) && override_thresholding)
        idx = 1;
    end
    if(isempty(idx))
        fidx(f) = fidx_orig(f);
    else
        chidx(f) = chansortidx(idx);
        fidx(f) = tmpidx(sortidx(1,chidx(f)));
    end
end

fidx_offset = fidx-fidx_orig;
%if(~all(fidx_offset == 0))
%	fprintf('freqtag_offset = [%d %d]\n',fidx_offset);
%end

if(nargout == 1)
	varargout = {fidx};
elseif(nargout == 2)
	varargout = {fidx, chidx};
elseif(nargout == 3)
    varargout = {fidx, chidx, fidx_offset};
end


