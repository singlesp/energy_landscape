function [event_raster] = event_related_events_raster(primary_event_idx, secondary_event_idx, windowsize)
%function [event_idx] = event_related_events(primary_event_idx, secondary_event_idx, windowsize)
%align secondary events around each primary event

[a b] = meshgrid(primary_event_idx(:), secondary_event_idx(:));
d = b-a;

event_raster = zeros(2*windowsize + 1, numel(primary_event_idx));
for i = 1:numel(primary_event_idx)
    event_raster(d(abs(d(:,i)) <= windowsize,i) + windowsize + 1,i) = 1;
end

