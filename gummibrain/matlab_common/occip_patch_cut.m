function P = occip_patch_cut(FV, cutvertidx, planevertidx)


planeverts = FV.vertices(planevertidx,:);

C = cross(planeverts(2,:)-planeverts(1,:),planeverts(3,:)-planeverts(1,:));
p = bsxfun(@minus,FV.vertices,planeverts(1,:))*(C.');
patchmask = (sign(p) == sign(mean(p(setdiff(cutvertidx,planevertidx)))));

A = mesh_adjacency_kj(FV.faces);
nei = vertex_neighbours(struct('faces',FV.faces));

pcut = shortest_path(cutvertidx,[],A,nei);
patchmask(pcut) = false;

patchfaces = FV.faces(all(patchmask(FV.faces),2),:);

ind = unique(reshape(patchfaces,[],1));
e = surfedge(patchfaces); %see iso2mesh toolbox

edgevert = unique(reshape(e,[],1));

P = struct('x',FV.vertices(ind,1),'y',FV.vertices(ind,2),'z',FV.vertices(ind,3), ...
    'ind', ind, 'vno', ind, 'edge', edgevert, 'npts', numel(ind));
