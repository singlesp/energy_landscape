function Vpad = padvolumeto(V,finalsize)

sz0=size(V);
%sz=max(sz0)*[1 1 1];
sz=finalsize;
if(numel(sz)==1)
    sz=sz*[1 1 1];
end

if(isequal(sz,sz0))
    Vpad=V;
    return
end

Vpad=zeros(sz(1),sz(2),sz(3));
Vpad(1:sz0(1),1:sz0(2),1:sz0(3))=V;
Vpad=circshift(Vpad,round((sz-sz0)/2));
