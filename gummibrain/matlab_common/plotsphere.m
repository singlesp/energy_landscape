function hsurf = plotsphere(pos,radius,cvalues,varargin)

if(nargin < 3 || isempty(cvalues))
    cvalues = [];
end

if(numel(varargin) > 1 && strcmpi(varargin{1},'resolution'))
    sphereres = varargin{2};
    varargin = {varargin{3:end}};
else
    sphereres = 20;
end

if(min(size(pos)) == 1)
	pos = reshape(pos,1,3);
end

[sx sy sz] = sphere(sphereres);
np = get(gca,'nextplot');
set(gca,'nextplot','add');
hsurf = zeros(size(pos,1),1);
for i = 1:size(pos,1)
    nr = min(i,numel(radius));
    hsurf(i) = surface(sx*radius(nr)+pos(i,1), sy*radius(nr)+pos(i,2), sz*radius(nr)+pos(i,3),varargin{:});
    if(~isempty(cvalues))
		nc = min(i,size(cvalues,1));
		if(size(cvalues,2) == 1)
			cv = repmat(cvalues(nc),size(get(hsurf(i),'CData')));
		elseif(size(cvalues,2) == 3)
			cv = repmat(reshape(cvalues(nc,:),[1 1 3]),size(get(hsurf(i),'CData')));
		end
        set(hsurf(i),'CData',cv);
    end
end
set(gca,'nextplot',np);