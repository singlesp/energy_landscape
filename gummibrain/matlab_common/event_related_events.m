function event_offset = event_related_events(primary_event_idx, secondary_event_idx, windowsize)
%function [event_idx] = event_related_events(primary_event_idx, secondary_event_idx, windowsize)
%align secondary events around each primary event


[a b] = meshgrid(primary_event_idx(:), secondary_event_idx(:));
d = b-a;
event_offset = d(:);

% 
% %figure; plot(d)
% if(nargin > 2 && ~isempty(windowsize))
%     event_offset = event_offset(abs(event_offset) <= windowsize);
% end

event_offset = {};
for i = 1:numel(primary_event_idx)
    event_offset{i} = d(d(:,i) >= 0 & d(:,i) < windowsize,i) + 1;
end