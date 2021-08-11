function [xyz faces] = cube(cubesize)
if(nargin < 1 || isempty(cubesize))
    cubesize = [1 1 1];
end

if(max(size(cubesize)) == 1)
    xsize = cubesize;
    ysize = cubesize;
    zsize = cubesize;
elseif(max(size(cubesize)) == 3)
    xsize = cubesize(1);
    ysize = cubesize(2);
    zsize = cubesize(3);
end

xyz = [0 0 0; 
    0 1 0;
    1 1 0;
    1 0 0;
    0 0 1; 
    0 1 1;
    1 1 1;
    1 0 1];

xyz = xyz*diag([xsize ysize zsize]);

% tri = [1 2 3;
%     1 3 4;
%     5 6 7;
%     5 7 8;
%     1 2 5;
%     2 6 5;
%     2 3 6;
%     3 7 6;
%     3 4 7;
%     4 8 7;
%     4 1 8;
%     1 5 8];

%[xyz tri] = meshcheckrepair(xyz,tri);
%used meshcheckrepair to make sure all faces going the right way
% tri = [1     8     5
%      4     8     1
%      4     7     8
%      3     7     4
%      3     6     7
%      2     6     3
%      2     5     6
%      1     5     2
%      5     8     7
%      5     7     6
%      3     4     1
%      2     3     1];
 
 % just do square patches
 faces = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8];

