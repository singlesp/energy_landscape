function v = generate_sphere_points(N,varargin)
% v = generate_sphere_points(N,[options...])
%
% Generates a set of N evenly-spaced points on a unit sphere 
%
% Options:
%   base: spiral(default), ico, oct
%       spiral: start with archimedes spiral of N points.  This is the only
%               method for generating exactly N points, and requires the
%               most repositioning iterations
%
%       ico:    start with icosahedron (12 points) and subdivide until
%               there are at least N points.  The following point counts
%               are possible: [12, 42, 162, 642, 2562, ...]=2+(12-2)*4^n
%
%       oct:    same as ico but start with octahedron (6 points)
%               Point counts: [6, 18, 66, 258, 1026, 4098, ...]=2+(6-2)*4^n
%
%
% uses external functions: 
% SpiralSampleSphere, mesh_subdivide,
% spm_matrix, spm_imatrix 

%default options
options=struct(...
    'base','spiral',...
    'maxiter',4000,...
    'dthresh',1e-10,...
    'fix_ends',true,...
    'verbose',false);

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
%parse options
input_opts=mergestruct(varargin{:});
fn=fieldnames(input_opts);
for f = 1:numel(fn)
    options.(fn{f})=input_opts.(fn{f});
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

switch(options.base)
    case 'spiral'
        [v, f]=SpiralSampleSphere(N);
    case 'oct'
        [v,f]=octahedron;
    case 'ico'
        [v,f]=icosahedron;
end

subN=0;
while size(v,1) < N
    subN=subN+1;
    [v,f]=mesh_subdivide(v,f);
    v=v./repmat(sqrt(sum(v.^2,2)),1,3);
end
Nnew=size(v,1);

if(subN>0)
    %make index increase in -z direction
    [~,ix]=sort(v(:,3),'descend');
    v=v(ix,:);
    f = unique(ix(f),'rows');
end

viewfix_idx=[];
if(options.fix_ends)
    viewfix_idx=[1 Nnew];
end

viewfix_mask=false(Nnew,1);
viewfix_mask(viewfix_idx)=true;

v0=v;
vtmp=v0;
dt=zeros(options.maxiter,1);
for i=1:options.maxiter
    G=point_repulsion(v);
    v=v+G;
    if(~isempty(viewfix_idx))
        v(viewfix_mask,:)=v0(viewfix_mask,:);
    end
    
    v=v./repmat(sqrt(sum(v.*v,2)),1,3);

    dt(i)=mean(sqrt(sum((v-vtmp).^2)));
    
    
    if(dt(i)<=options.dthresh)
        if(options.verbose)
            fprintf('delta<=%.4d after %d iterations\n',options.dthresh,i);
        end
        break
    end
    vtmp=v;

end
if(options.verbose && options.maxiter>0 && dt(end)>options.dthresh)
    fprintf('delta=%.4d after %d iterations\n',dt(end),options.maxiter);
end

if(~options.fix_ends)
    A=([v ones(Nnew,1)]\[v0 ones(Nnew,1)])';
    p=spm_imatrix(A);
    p=[0 0 0 p(4:6) 1 1 1];
    A=spm_matrix(p);

    v = affine_transform(A,v);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function vgrad = point_repulsion(vxyz)
N=size(vxyz,1);
x=repmat(vxyz(:,1),1,N);
y=repmat(vxyz(:,2),1,N);
z=repmat(vxyz(:,3),1,N);

dx=x-x';
dy=y-y';
dz=z-z';

d2=dx.*dx+dy.*dy+dz.*dz;
d2(1:(N+1):end)=1; %set diagonals to 1 to avoid div. by 0
vgrad=[sum(dx./d2,2) sum(dy./d2,2) sum(dz./d2,2)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,f] = icosahedron
t = (1+sqrt(5)) / 2;
% create vertices
v = [-1, t, 0; % v1
      1, t, 0; % v2
     -1,-t, 0; % v3
      1,-t, 0; % v4
      0,-1, t; % v5
      0, 1, t; % v6
      0,-1,-t; % v7
      0, 1,-t; % v8
      t, 0,-1; % v9
      t, 0, 1; % v10
     -t, 0,-1; % v11
     -t, 0, 1];% v12
% normalise vertices to unit size
v = v./norm(v(1,:));

% create faces
f = [ 1,12, 6; % f1
      1, 6, 2; % f2
      1, 2, 8; % f3
      1, 8,11; % f4
      1,11,12; % f5
      2, 6,10; % f6
      6,12, 5; % f7
     12,11, 3; % f8
     11, 8, 7; % f9
      8, 2, 9; % f10
      4,10, 5; % f11
      4, 5, 3; % f12
      4, 3, 7; % f13
      4, 7, 9; % f14
      4, 9,10; % f15
      5,10, 6; % f16
      3, 5,12; % f17
      7, 3,11; % f18
      9, 7, 8; % f19
     10, 9, 2];% f20


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v,f] = octahedron
v=[0 0 1
    1 0 0
    0 1 0
    -1 0 0
    0 -1 0
    0 0 -1];
    
f=[1 3 2
    1 4 3
    1 5 4
    1 2 5
    6 2 3
    6 3 4
    6 4 5
    6 5 2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function v = affine_transform(M, v)
%v = affine_transform(M, v)
% M: 4x4
% input: Nx3
% output: Nx3

v = M * [v.'; ones(1, size(v, 1))];
v = v(1:3,:)';
