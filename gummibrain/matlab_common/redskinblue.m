function C = redskinblue(m)
%     n1  n2   m/2
%r = 0->0->.25->.5    .5->1->.5
%g = 0->0->.25->.5     .5->0->0
%b = .5->1->.75->.5   .5->0->0

skincol = .75*[0.9373 0.8157 0.8118];
p = [0 .2 .45 .55 .8 1];
c = [0 0 .6; 0 0 1; skincol; skincol; 1 0 0; .5 0 0];
C = GenerateColormap(p,c,m);