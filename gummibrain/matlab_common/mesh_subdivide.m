function [v,f] = mesh_subdivide(v,f,n)
%
% Adapted from icosphere.m by:
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK

if(~exist('n','var') || isempty(n))
    n=1;
end

% recursively subdivide triangle faces
for gen = 1:n
    f_ = zeros(size(f,1)*4,3);
    for i = 1:size(f,1) % for each triangle
        tri = f(i,:);
        % calculate mid points (add new points to v)
        [a,v] = getMidPoint(tri(1),tri(2),v);
        [b,v] = getMidPoint(tri(2),tri(3),v);
        [c,v] = getMidPoint(tri(3),tri(1),v);
        % generate new subdivision triangles
        nfc = [tri(1),a,c;
               tri(2),b,a;
               tri(3),c,b;
                    a,b,c];
        % replace triangle with subdivision
        idx = 4*(i-1)+1:4*i;
        f_(idx,:) = nfc;
    end
    f = f_; % update 
end

% % remove duplicate vertices
[v,b,ix] = unique(v,'rows'); clear b % b dummy / compatibility
% % reassign faces to trimmed vertex list and remove any duplicate faces
f = unique(ix(f),'rows');

end

function [i,v] = getMidPoint(t1,t2,v)
%GETMIDPOINT calculates point between two vertices
%   Calculate new vertex in sub-division and normalise to unit length
%   then find or add it to v and return index
%
%   Wil O.C. Ward 19/03/2015
%   University of Nottingham, UK

% get vertice positions
p1 = v(t1,:); p2 = v(t2,:);
% calculate mid point (on unit sphere)
pm = (p1 + p2) ./ 2;

% add to vertices list, return index
i = size(v,1) + 1;
v = [v;pm];

end
