function gfp = GlobalFieldPower(eeg)
% gfp = GlobalFieldPower(eeg)
%
% eeg = ExT matrix where each column is a topography
% returns global field power for each topography

[nr nc] = size(eeg);
if(nr==1 || nc==1)
    eeg = eeg(:);
end

eeg = eeg';

[num_timepoints num_elec] = size(eeg);

gfp = zeros(num_timepoints,1);

for t = 1:num_timepoints
    dsq = distance(eeg(t,:),eeg(t,:)).^2;
    gfp(t) = sum(dsq(:));
end

gfp = sqrt(gfp/(2*num_elec));
%gfp = sqrt(gfp/(num_elec*num_elec));