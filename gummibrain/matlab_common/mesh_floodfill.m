function meshmask = mesh_floodfill(meshvals_in, conn, startidx, maxiter,bequiet)

if(~exist('maxiter','var') || isempty(maxiter))
    maxiter = 1000;
end

if(~exist('bequiet','var') || isempty(bequiet))
    bequiet = false;
end

meshmask = ones(size(meshvals_in));
val = meshvals_in(startidx);
meshmask(meshvals_in ~= val) = 0;

if(~isempty(which('mesh_floodfill_double')))
    meshmask_out = mesh_floodfill_double(+meshmask,conn,startidx,maxiter,bequiet);
    meshmask = (meshmask_out == 2);
else
    %%%% if mesh_floodfill_double not compiled
    fprintf('mesh_floodfill_double.c not compiled.  Using slower Matlab code\n');
    idx = startidx;

    meshmask(idx) = 2;
    meshcount = 1;
    for i = 1:maxiter
        newidx = cat(1,conn{idx});
        newidx = newidx(meshmask(newidx) == 1);
        meshmask(newidx) = 2;
        idx = newidx;
        tmp_meshcount = sum(meshmask == 2);
        if(tmp_meshcount == meshcount)
            break;
        end
        meshcount = tmp_meshcount;
        %fprintf('%d\n',meshcount);
    end
    if(~bequiet)
        if(i == maxiter)
            fprintf('Max iterations (%d) exceeded.  Quitting.\n\n',maxiter);
        else
            fprintf('Filled %d vertices in %d iterations.\n\n',meschount,i);
        end
    end
    meshmask = (meshmask == 2);
end