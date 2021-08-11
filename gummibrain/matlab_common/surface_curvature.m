function curv = surface_curvature(faces,verts)

norms = surfnormals(verts,faces);

num_verts = size(verts,1);

conn = meshconn(faces,max(faces(:)));

curv = zeros(num_verts,1);

for i = 1:num_verts
    nei = [conn{i}];

    %if(numel(nei) == 0)
    %    continue;
    %end
    
    nei_dots = norms(nei,:) * norms(i,:);
    curv(i) = mean(nei_dots);
end
% 
% for i = 1:num_verts
%     nei = faces(any(faces == i,2),:);
%     nei = unique(nei(nei ~= i));
%     
%     if(numel(nei) == 0)
%         continue;
%     end
%     
%     nei_norms = norms(nei,:);
%     nei_dots = nei_norms * norms(i,:)';
% 
%     curv(i) = mean(nei_dots);
% end