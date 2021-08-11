function HideFromLegend(h)

set(get(get(h,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude line from legend
