
function [bad_elec good_elec] = GetBadElectrodesNS64(subj)


num_elec = 68;
bad_elec = [];
bad_elec = [33 43 65:68];
%bad_elec = [31 32]; %EOG, ECG
%bad_elec = [10 11 84 85 108 110 111 120 129:132]; %missing elec (12)

switch subj
    case 'ly20110407'
        bad_elec = [bad_elec];
    case 'allsubj'
    case 'allsubj_interp'
end

good_elec = 1:num_elec;
good_elec(bad_elec) = [];