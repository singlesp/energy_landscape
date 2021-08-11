function cmap = colormap_retinotopy_wedge(n)
if(~exist('n','var') || isempty(n))
	n=128;
end

p = 0:.125:1;

cw = [0 0 1;
	0 .75 .75; 
	0 .5 0; 
	.5 .75 0;
	1 1 0;
	1 .5 0;
	1 0 0;
	1 0 1;
	0 0 1];

cmap = GenerateColormap(p,cw,n);