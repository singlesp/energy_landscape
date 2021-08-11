function nei = gridneighbors(gridsize)

idx = 1:prod(gridsize);
[i j k] = ind2sub(gridsize,idx);
