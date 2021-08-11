function h = latexticks(ax)

if(nargin == 0 || isempty(ax))
    ax=gca;
end

%%% Get tick mark positions
yTicks = get(ax,'ytick');
yTickLabel = cellstr(get(ax,'yticklabel'));
xTicks = get(ax,'xtick');
xTickLabel = cellstr(get(ax,'xticklabel'));

set(ax,'yticklabel',[],'xticklabel',[]); %Remove tick labels
axlim = axis; %Get left most x-position
axpos = get(ax,'position');
axpos2 = getpixelposition(ax);

%fontsize = 14;
fontsize = get(ax,'fontsize');


h = [];

fpp = axpos(3:4)./axpos2(3:4);
HorizontalOffset = .5*fpp(1);
VerticalOffset = 30*fpp(2);

%%% Reset the ytick labels in desired font
for i = 1:length(yTicks)
    %Create text box and set appropriate properties
    h(end+1) = text(axlim(1) - HorizontalOffset,yTicks(i),[yTickLabel{i}],...
        'HorizontalAlignment','Right','fontsize',fontsize);
end

for i = 1:length(xTicks)
    %Create text box and set appropriate properties
    
    if(strcmpi(get(ax,'YDir'),'reverse'))
        y = axlim(4);
    else
        y = axlim(3);
    end
    y = y-VerticalOffset;
    h(end+1) = text(xTicks(i),y - VerticalOffset, [xTickLabel{i}],...
        'HorizontalAlignment','Center','VerticalAlignment','Top','fontsize',fontsize);
end
%set(ax,'fontsize',fontsize);
