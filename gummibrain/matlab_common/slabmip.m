function img = slabmip(vol,ax,slice,slabthick,minmax)

if(~exist('minmax','var') || isempty(minmax))
    minmax='min';
end

slabfun=@(x)min(x,[],3);
switch lower(minmax)
    case 'max'
        slabfun=@(x)max(x,[],3);
    case 'mean'
        slabfun=@(x)mean(x,3);
    case 'min'
    otherwise
end
sz=size(vol);
vol=permute(vol,[setdiff(1:3,ax) ax]);
img=slabfun(vol(:,:,max(1,min(sz(ax),slice+(1:slabthick)))));
