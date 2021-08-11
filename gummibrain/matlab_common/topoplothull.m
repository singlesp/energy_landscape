function varargout = topoplothull(varargin)

topoplot(varargin{:},'style','both','electrodes','labelpoint');
ax = gca;

hlabels = findobj(ax,'type','text');
labelpos = cell2mat(get(hlabels,'position'));
labelpos = labelpos(:,1:2);

xl = get(gca,'xlim');
yl = get(gca,'ylim');
limpos = [xl([1 2 2 1]); yl([1 1 2 2])]';

vert = [labelpos; limpos];
tri = delaunay(vert);
tri = tri(~all(tri <= size(labelpos,1),2),:);

ax2 = copyobj(gca,gcf);
set(gcf,'currentaxes',ax2);
delete(ax);
delete(get(ax2,'children'));

hold on;
varargout = topoplot(varargin{:});
h=patch('vertices',vert,'faces',tri,'linestyle','none','facecolor','w','facealpha',1);