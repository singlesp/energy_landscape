function densevals = transfer2densesurf(denseverts,densefaces,vals,interptype)
% densevals = transfer2densesurf(verts,denseverts,densefaces,vals,interptype)
%
% For each new vertex in dense surface, identify the TWO original vertices
%  from which it was generated, compute the distance to each of those
%  original vertices, and use either the value from the nearest of those
%  vertices (interptype='nearest') or a distance-weighted average
%  (interptype='linear')
%
% This function assumes the dense surface was generated the way
%  freesurfer mris_mesh_subdivide is used: New vertices are halfway along
%  each edge of original surface and denseverts is [origverts; newverts]
%
% For repeat operations, you can generate a "lookup" as follows:
% >> denselookup=transfer2dense(verts,denseverts,densefaces,size(verts,1));
% >> densevals=vals(denselookup) % if vals is an Nx1 vector, or
% >> densevals=vals(denselookup,:) % if vals is a NxT matrix
%
% Inputs:
%   denseverts:     Dx3 vertex positions from dense surface
%   densefaces:     Fx3 faces from dense surface
%   vals:           NxT vertex values on original surface
%          or       1x1 scalar specifying the number of vertex on original
%                     surface. If vals=<numorig>, use 1:<numorig> to 
%                     create lookup (in which case interptype is changed 
%                     to 'nearest')
%   interptype:     'linear' or 'nearest' (default)
%
% Outputs:
%   densevals:      MxT vertex values on dense surface
%
if(~exist('interptype','var') || isempty(interptype))
    interptype='nearest';
end

if(numel(vals)==1)
    oldN=vals;
    vals=(1:oldN)';
    interptype='nearest';
elseif(numel(vals)>1)
    if(size(vals,1)==1)
        vals=vals.';
    end
    oldN=size(vals,1);
else
    error('input vals cannot be empty');
end


edgesd=[densefaces(:,[1 2]); densefaces(:,[1 3]); densefaces(:,[2 3])];

newidx=oldN+1:size(denseverts,1);
vnew=denseverts(newidx,:);

minedge=min(edgesd,[],2);
maxedge=max(edgesd,[],2);

newedge=edgesd(minedge<=oldN & maxedge>oldN,:);

[~,idx]=sort(max(newedge,[],2));
newedge=newedge(idx,:);
origidx=min(newedge,[],2);

origidx=reshape(origidx,4,[])';
orig1=min(origidx,[],2);
orig2=max(origidx,[],2);

d1=sum((vnew-denseverts(orig1,:)).^2,2);
d2=sum((vnew-denseverts(orig2,:)).^2,2);

if(isequal(lower(interptype),'nearest'))
    densevals=vals(orig1,:);
    densevals(d2<d1,:)=vals(orig2(d2<d1),:);
else
    d1=sqrt(d1);
    d2=sqrt(d2);
    d1frac=d1./(d1+d2);
    densevals=vals(orig1,:).*(1-d1frac)+vals(orig2,:).*(d1frac);
end
densevals=[vals; densevals];
