function goodfaces = RemoveLongEdges(faces,verts, maxlength)

nf = size(faces,1);

edge_faces = repmat((1:nf)',3,1);

edgelen = edgelengths(faces,verts);

bad_edges = find(edgelen > maxlength);
bad_faces_idx = unique(edge_faces(bad_edges));
goodface_bool = true(nf,1);
goodface_bool(bad_faces_idx) = false;

goodfaces = faces(goodface_bool,:);