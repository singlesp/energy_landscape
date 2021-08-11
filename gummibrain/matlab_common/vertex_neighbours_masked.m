function Ne = vertex_neighbours_masked(faces,vertmask)

maskedfaces=faces(any(vertmask(faces),2),:);

maskedverts_and_neighbors=false(size(vertmask));
maskedverts_and_neighbors(maskedfaces)=1;

idxmasked=find(maskedverts_and_neighbors);
ordermasked=(1:numel(idxmasked))';
vtmp=zeros(size(maskedverts_and_neighbors));
vtmp(idxmasked)=ordermasked;
ftmp=vtmp(maskedfaces);

Ne=vertex_neighbours(struct('faces',ftmp,'vertices',zeros(size(idxmasked,1),3)));
Ne=Ne(vtmp(vertmask,:));

Ne=cellfun(@(i)idxmasked(i),Ne,'uniformoutput',false);
