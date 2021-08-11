function cmap = colormap_retinotopy_wedge2(n)
if(~exist('n','var') || isempty(n))
	n=128;
end

%p = [0 1/6 1/3 3/6 4/6 5/6 1];
p=linspace(0,1,7);
cw = [0 0  1
 0 0 .5
 0 1 0
 .5 0 0
 1 0 0
 0 0 0
 0 0 1];

cmap = GenerateColormap(p,cw,n);