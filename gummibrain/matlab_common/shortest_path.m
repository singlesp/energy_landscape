function p = shortest_path(v1,v2,faces_or_adjacency,nei)
if(issparse(faces_or_adjacency) || size(faces_or_adjacency,1) == size(faces_or_adjacency,2))
    A = faces_or_adjacency;
else
    A = mesh_adjacency_kj(faces_or_adjacency);
end


if(numel(v1) > 1)
    p = v1(1);
    for n = 1:numel(v1)-1
        ptmp = shortest_path(v1(n),v1(n+1),A,nei);
        p = [p ptmp(2:end)];
    end
    return;
end

G = dijkstra_kj(A,v2);
if(G(v1) == 0)
    p = [];
    return;
end

p = zeros(1,G(v1)+1);
cur = v1;
p(1) = cur;
for i = 1:G(v1)
    nt = nei{cur};
    [~,idx] = min(G(nt));
    cur = nt(idx);
    %nt = find(A(cur,:) > 0);
    %[m,idx] = min(G(nt));
    %cur = nt(idx);
    G(cur) = inf;
    p(i+1) = cur;
end
