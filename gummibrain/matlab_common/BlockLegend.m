function h = BlockLegend(fig,names,colors,varargin)

hf=figure;
for i = 1:numel(names)
    bar(i,1,'facecolor',colors(i,:));
    hold on;
end

hl=legend(names,varargin{:});
h=copyobj(hl,fig);
close(hf);

locidx = StringIndex(varargin,'location');
if(isempty(locidx) || numel(varargin) == locidx)
    return;
end

locstr = lower(varargin{locidx+1});

lp=get(h,'position');
lw = lp(3);
lh = lp(4);

x = 0; %southwest
y = 0;

if(strfind(locstr,'west'))
    x = 0;
end

if(strfind(locstr,'east'))
    x = 1-lw;
end

if(strfind(locstr,'north'))
    y = 1-lh;
end

if(strfind(locstr,'south'))
    y = 0;
end

set(h,'position', [x y lw lh]);
