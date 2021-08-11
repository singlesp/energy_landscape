function ax = gaa(fig)
if(nargin < 1 || isempty(fig))
    fig = gcf;
end
ax = findobj(fig,'type','axes');
axtag = get(ax,'tag');
axcb = strcmpi(axtag,'colorbar');
axleg = strcmpi(axtag,'legend');
axtitle = strcmpi(axtag,'suptitle');
ax = ax(~axcb & ~axleg & ~axtitle);