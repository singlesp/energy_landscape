function axisloosen(ax,amt)

if(~exist('ax','var') || isempty(ax))
    ax = gca;
end

if(~exist('amt','var') || isempty(amt))
    amt = 1.1;
end

for i = 1:numel(ax)
    xl = get(ax(i),'xlim');
    yl = get(ax(i),'ylim');
    xl = amt*(xl-mean(xl))+mean(xl);
    yl = amt*(yl-mean(yl))+mean(yl);
    set(ax(i),'xlim',xl);
    set(ax(i),'ylim',yl);
end