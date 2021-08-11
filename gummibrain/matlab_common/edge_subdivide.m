function [v,f] = edge_subdivide(v,f,n)

vnew=v;
if(size(f,2)==1)
    edges=[f(1:end-1) f(2:end)];
elseif(size(f,2)==2)
    edges=f;
elseif(size(f,2)==2)
    edges=zeros(size(f,1)*3,2);
    edges(1:3:end,:)=f(:,[1 2]);
    edges(2:3:end,:)=f(:,[2 3]);
    edges(3:3:end,:)=f(:,[1 3]);
else
    error('Cannot handle faces with more than 3 vertices');
end

for i = 1:n
    vnew=[vnew; (v(edges(:,1),:)+v(edges(:,2),:))/2];
end
