function [lap_vals lap_surf_idx] = sphere_laplacian(elec_vals, elec_pos, surf_vals,lap_shape, lap_radius)

if(nargin < 3)
    surf_vals = [];
end

if(nargin < 4)
    lap_shape = [];
end

if(nargin < 5)
    lap_radius = [];
end

sph_file = 'sphere_surf_simp.mat';
Sph = load(sph_file);

[lap_vals lap_idx] = surface_laplacian(elec_vals, elec_pos, surf_vals, Sph.POS, Sph.NORMS, lap_shape, lap_radius);

if(nargout > 1)
    lap_surf_idx = lap_idx;
end
