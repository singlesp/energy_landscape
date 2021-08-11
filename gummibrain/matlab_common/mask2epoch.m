function [epochs, epoch_mat] = mask2epoch(mask)

newepochtimes = [1; find(diff(mask(:)) ~= 0) + 1];
epoch_times = [mask(newepochtimes) newepochtimes];
epoch_times = sortrows(epoch_times,2);

if(isempty(epoch_times))
    epochs = [];
else
    epoch_times = [epoch_times [diff([epoch_times(:,2); size(mask,1)+1])]];
    epochs = struct('type',num2cell(epoch_times(:,1)),...
        'latency',num2cell(epoch_times(:,2)),...
        'duration',num2cell(epoch_times(:,3)));
end
if(nargout > 1)
    epoch_mat = epoch_times;
end
