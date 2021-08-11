function [vv,ff] = octahedron()
vv=[0 0 1
    1 0 0
    0 1 0
    -1 0 0
    0 -1 0
    0 0 -1];
    
ff=[1 3 2
    1 4 3
    1 5 4
    1 2 5
    6 2 3
    6 3 4
    6 4 5
    6 5 2];
