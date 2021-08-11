function [M_laplace lap_center_idx lap_surround_idx lap_shape lap_radius M_evalpos] = sphere_laplacian_transfer(splinefile, lap_shape, lap_radius, eval_pos)
%%%%%%%%%%%% build surface laplacian transfer matrix for this specific electrode arrangement

if(nargin < 2)
    lap_shape = [];
end

if(nargin < 3)
    lap_radius = [];
end

if(nargin < 4)
    eval_pos = [];
end

sph_file = 'sphere_surf_simp.mat';
Sph = load(sph_file);

Sp = load(splinefile,'newElect');
elec_pos = Sp.newElect;

num_eval = size(eval_pos,1);

[c r] = spherefit(elec_pos);
elec_pos = [elec_pos; eval_pos];
elec_pos = elec_pos - repmat(c',size(elec_pos,1),1);
elec_pos = normrows(elec_pos);

eval_pos = elec_pos(end-num_eval+1:end,:);
elec_pos = elec_pos(1:end-num_eval,:);

[M_laplace lap_center_idx lap_surround_idx lap_shape lap_radius M_evalpos] = surface_laplacian_transfer(splinefile, elec_pos, Sph.POS, Sph.NORMS, lap_shape, lap_radius, eval_pos);
