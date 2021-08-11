function EEG = icasort(EEG)

nchans = size(EEG.icaweights,2);
winv_strength = sqrt(sum(EEG.icawinv.^2,1)/nchans);
w_strength = sqrt(sum(EEG.icaweights.^2,2)'/nchans);
ica_strength = winv_strength.*w_strength;

EEG.icawinv = EEG.icawinv*diag(1./winv_strength);
EEG.icaweights = diag(ica_strength./w_strength)*EEG.icaweights;
EEG.data = double(EEG.data);
EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data(EEG.icachansind,:);

%
% compute variances without backprojecting to save time and memory -sm 7/05
%

%meanvar = sqrt(sum(EEG.icawinv.^2).*sum((EEG.icaact').^2)./(numel(EEG.icachansind)*size(EEG.data,2)-1)); % from Rey Ramirez 8/07
%meanvar = var(EEG.icaact,[],2);
%
%%%%%%%%%%%%%% Sort components by mean variance %%%%%%%%%%%%%%%%%%%%%%%%
%
[~, windex] = sort(ica_strength,'descend');

EEG.icaact = EEG.icaact(windex,:);
EEG.icawinv = EEG.icawinv(:,windex);
EEG.icaweights = EEG.icaweights(windex,:);

