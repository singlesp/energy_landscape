function eloc = electrode_dat2eeglab(datfile,missing_elec_names,recenter,remove_landmarks,alt_landmarks)
% eloc = electrode_dat2eeglab(datfile,missing_elec_names)

if(nargin < 2)
    missing_elec_names = [];
end

if(nargin < 3)
    recenter = true;
end

if(nargin < 4)
    remove_landmarks = true;
end

if(nargin < 5)
    alt_landmarks = [];
end

missing_elec_idx = [];
if(isnumeric(missing_elec_names))
    missing_elec_idx = missing_elec_names;
    missing_elec_names = {'veo','heo','ekg','emg'};
end

M = importdata(datfile);

elec_name = M.textdata;

bad_name = {'centroid','ref.'};
if(remove_landmarks)
    bad_name = {bad_name{:}, 'nasion','left','right'};
elseif(~isempty(alt_landmarks))
    M.data(1:3,:) = M.data(3+alt_landmarks,:);
end

elec_name_sp = strcat('_',elec_name,'_'); %pad with nonsense character so avoid partial matches (11 matching 115)
badcat = strcat('_',bad_name,'_|');
badcat = [badcat{:}];
good_elec = cellfun(@isempty,regexpi(elec_name_sp,badcat));
%good_elec = find(cellfun(@isempty,regexpi(elec_name_sp,badcat)));
%bad_elec = find(cellfun(@numel,regexpi(elec_name_sp,badcat)));

elec_name = {elec_name{good_elec}};
data_xyz = M.data(good_elec,2:4);

if(~isempty(missing_elec_names))
    bad_name = missing_elec_names;
    elec_name_sp = strcat('_',elec_name,'_'); %pad with nonsense character so avoid partial matches (11 matching 115)
    badcat = strcat('_',bad_name,'_|');
    badcat = [badcat{:}];
    good_elec = cellfun(@isempty,regexpi(elec_name_sp,badcat)); %overwrite
    %good_elec = find(cellfun(@isempty,regexpi(elec_name_sp,badcat)));
    %bad_elec = find(cellfun(@numel,regexpi(elec_name_sp,badcat)));
end

if(~isempty(missing_elec_idx))
    good_elec(missing_elec_idx) = false;
end

data_xyz = data_xyz(good_elec,:);
elec_name = {elec_name{good_elec}};


if(remove_landmarks)
    [c r] = spherefit(data_xyz);
else
    [c r] = spherefit(data_xyz(4:end,:));
end

if(recenter)
    data_xyz = data_xyz - repmat(c',size(data_xyz,1),1);
end
%data_xyz = data_xyz ./ r;

eloc = struct('Y',num2cell(-data_xyz(:,1)),...
    'X',num2cell(data_xyz(:,2)),...
    'Z',num2cell(data_xyz(:,3)),...
    'labels',elec_name');

eloc = convertlocs(eloc, 'cart2all');