function bbox = getbbox(X,pad)

sz=size(X);
if(~exist('pad','var') || isempty(pad))
    pad=zeros(1,numel(sz));
end

bbox=zeros(2,numel(sz));

for d=1:numel(sz)
    Xtmp=any(reshape(permute(X,[setdiff(1:3,d) d]),[],sz(d)),1);
    bbox(1,d)=find(Xtmp,1,'first');
    bbox(2,d)=find(Xtmp,1,'last');
    bbox(:,d)=bbox(:,d)+pad(d)*[-1; 1];
    bbox(:,d)=min(sz(d),max(1,bbox(:,d)));
end
%bbox=bbox+[-pad; pad];
