function V_diffuse = mesh_diffuse(V,Ne,iter)
%V_diffuse = mesh_diffuse(V,Ne,iter)
%Use Ne=vertex_neighbours(struct('faces',F,'vertices','V'));

if(numel(V) < max(cat(1,Ne{:})))
    error('not enough vertices for this mesh neighbor list');
end

if(~isreal(V))
    %Vr = mesh_diffuse(real(V),Ne,iter);
    %Vi = mesh_diffuse(imag(V),Ne,iter);
    %V_diffuse = complex(Vr,Vi);
    
    Vri = mesh_diffuse([real(V) imag(V)],Ne,iter);
    V_diffuse = complex(Vri(:,1:end/2),Vri(:,end/2+1:end));

    return;
end

V_diffuse = V;

if(size(V,1) < size(V,2))
    V_diffuse = V_diffuse.';
end

V_diffuse = mesh_diffuse_double(V_diffuse,Ne,iter);
return;
%%

for i = 1:iter
   for j = 1:size(V_diffuse,2)
       V_diffuse(:,j) = mesh_diffuse_double(V_diffuse(:,j),Ne);
   end
end
