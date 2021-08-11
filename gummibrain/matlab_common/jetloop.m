function cmap = jetloop(n)

p = [0 .125 .25 .50 .75 1 1.125];
c = [0 0 .5; 0 0 1; 0 1 1; 1 1 0; 1 0 0; .5 0 0; 0 0 .5];

cmap = GenerateColormap(p,c,n);
