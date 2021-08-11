function cmap = colormap_retinotopy_ring(n)
if(~exist('n','var') || isempty(n))
	n=128;
end

p = 0:.125:1;

cr = [0 .5 0;
	0 .75 .75;
	0 0 1;
	1 0 1;
	1 0 0;
	1 .5 0;
	1 1 0;
	.5 .75 0;
	0 .5 0];

cmap = GenerateColormap(p,cr,1000);
cmap = circshift(cmap,-round(.025*size(cmap,1)));
cmap = GenerateColormap(linspace(0,1,size(cmap,1)),cmap,n);