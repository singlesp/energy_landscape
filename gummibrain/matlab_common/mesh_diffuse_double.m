function V_diffuse = mesh_diffuse_double(V,Ne,iter)

V_diffuse = zeros(size(V));
for j=1:iter
    for i = 1:numel(Ne)
        V_diffuse(i) = mean(V(Ne{i}));
        V=V_diffuse;
    end
end
