function A = mesh_adjacency_kj(faces,numverts)
% From: triangulation2adjacency
%   Copyright (c) 2005 Gabriel Peyr

if(~exist('numverts','var') || isempty(numverts))
    numverts=max(faces(:));
end

f = double(faces);

A = sparse([f(:,1); f(:,1); f(:,2); f(:,2); f(:,3); f(:,3)], ...
           [f(:,2); f(:,3); f(:,1); f(:,3); f(:,1); f(:,2)], ...
           1.0);
       
% avoid double links
A = double(A>0);

if(size(A,1)<numverts || size(A,2)<numverts)
    Atmp=A;
    A=sparse(numverts,numverts);
    A(1:size(Atmp,1),1:size(Atmp,2))=Atmp;
end
