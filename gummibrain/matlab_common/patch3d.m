function h = patch3d(tri,xyz,varargin)

h = trisurf(tri,xyz(:,1),xyz(:,2),xyz(:,3),varargin{:});