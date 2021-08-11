function verts_new = mesh_smooth_vertices(verts, conn, weights, iter, alpha, beta)

if(~exist('iter','var') || isempty(iter))
    iter = 1;
end
if(~exist('alpha','var') || isempty(alpha))
    alpha = 1;
end
if(~exist('beta','var') || isempty(beta))
    beta = .5;
end
if(~exist('weights','var') || isempty(weights))
    weights = ones(size(verts,1),1);
end
verts_new = mesh_smooth_double(verts,conn,weights,iter,alpha,beta);