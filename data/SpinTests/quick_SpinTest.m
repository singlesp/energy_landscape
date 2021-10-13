function[rho_emp, p_perm] = quick_SpinTest(x, y, atlas_name, names, corr_type)
%% NB must be cortical only!
%Exception is Lausanne

if not(size(x,1) == size(y,1))
    disp('ERROR: mismatched input sizes!')
else
    num_rois = size(x,1);
end

if not(exist('names', 'var'))
    names = {'X', 'Y'};
end

if not(exist('corr_type', 'var'))
    corr_type = 'Spearman';
end

if strcmp(atlas_name, 'Lausanne')
    load(['/Users/sps253/Documents/ROI_maps/Lausanne_', num2str(num_rois), '_subnetworks.mat'])
    cortex = find(subnetworks<=7);
    x = x(cortex);
    y = y(cortex);
end
    
%Load pre-computed rotated maps
load(['/Users/sps253/Documents/MATLAB/SpinTests/rotated_maps/rotated_', atlas_name, '_', num2str(num_rois), '.mat'])

if min(min(perm_id)) == 0
    perm_id = perm_id+1; %add 1 since starting from 0
end

 [p_perm, rho_emp, rho_null_xy, rho_null_yx] = perm_sphere_p_al857(x,y,perm_id);

%Find the top and bottom 5%
all_rhos = [rho_null_xy; rho_null_yx'];
sorted_rhos = sort(all_rhos, 'descend');
top_cutoff = sorted_rhos(numel(sorted_rhos).*0.05);
bottom_cutoff = sorted_rhos(numel(sorted_rhos).*0.95);

%Plot histogram with cutoffs
histo_handle = figure; histogram(all_rhos, 40, 'Normalization', 'probability')
hold on
yax=ylim; plot([rho_emp, rho_emp], [yax(1), yax(2)], 'r', 'LineWidth', 2)
yax=ylim; plot([top_cutoff, top_cutoff], [yax(1), yax(2)], 'b', 'LineWidth', 2)
yax=ylim; plot([bottom_cutoff, bottom_cutoff], [yax(1), yax(2)], 'b', 'LineWidth', 2)

ylabel('Probability of occurrence')
xlabel([corr_type, ' correlation'])
title({[names{1}, ' vs ', names{2}];  ['r = ', num2str(rho_emp), ' p-perm = ', num2str(p_perm)]})


        