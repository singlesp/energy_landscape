function nextfig(do_hidden)
% cycle through figure windows
% nextfig
% nextfig(true) ==> include non-handle windows


if ~exist('do_hidden','var') || isempty(do_hidden)
    do_hidden=false;
end

if(do_hidden)
    figs=findall(0,'type','figure');
else
    figs=findobj(0,'type','figure');
end
if(isempty(figs))
    return;
end

figs=figs(end:-1:1)';
curfig=gcf;
curfigidx=find(figs==curfig,1,'first');
if(isempty(curfigidx))
    curfigidx=0;
end


nextfigidx=mod(curfigidx,numel(figs))+1;

figure(figs(nextfigidx));