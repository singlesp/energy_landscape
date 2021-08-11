function [M_laplace lap_center_idx lap_surround_idx lap_shape lap_radius M_evalpos] = surface_laplacian_transfer(splinefile, elec_pos, surf_verts, surf_norms, lap_shape, lap_radius, eval_pos)
%%%%%%%%%%%% build surface laplacian transfer matrix for this specific electrode arrangement

if(nargin < 5)
    lap_shape = [];
end

if(nargin < 6)
    lap_radius = [];
end

if(nargin < 7)
    eval_pos = [];
end

Sp = load(splinefile);

[~, Smat] = EvaluateSplineSurface([],Sp.G,Sp.gx);

if(isempty(elec_pos))
    %elec_pos = [[Sp.Xe] [Sp.Ye] [Sp.Ze]];
    elec_pos = Sp.newElect;
end

num_eval = size(eval_pos,1);

[~, lap_center_idx lap_surround_idx lap_shape lap_radius] = surface_laplacian([], elec_pos, [], surf_verts, surf_norms, lap_shape, lap_radius,eval_pos);

lapeval_center_idx = lap_center_idx(end-num_eval+1:end);
lap_center_idx = lap_center_idx(1:end-num_eval);

lapeval_surround_idx = lap_surround_idx(:,end-num_eval+1:end);
lap_surround_idx = lap_surround_idx(:,1:end-num_eval);
% 
% figure;
% plot3(surf_verts(:,1),surf_verts(:,2),surf_verts(:,3),'b.','markersize',1);
% hold on;
% plot3(surf_verts(lapeval_center_idx,1),surf_verts(lapeval_center_idx,2),surf_verts(lapeval_center_idx,3),'b.');
% plot3(surf_verts(lapeval_surround_idx,1),surf_verts(lapeval_surround_idx,2),surf_verts(lapeval_surround_idx,3),'r.');
% 
% plot3(surf_verts(lap_center_idx,1),surf_verts(lap_center_idx,2),surf_verts(lap_center_idx,3),'k.');

num_elec = numel(Sp.Xe);
M_laplace = zeros(num_elec,num_elec);

for n = 1:num_elec
    vals = zeros(num_elec,1);
    vals(n) = 1;
    
    Svals = EvaluateSplineSurface(vals,Smat);
    lap_vals = vals-mean(Svals(lap_surround_idx),1)';
    
    M_laplace(:,n) = lap_vals;
end

if(~isempty(eval_pos))
    M_evalpos = zeros(num_eval,num_elec);

    for n = 1:num_elec
        vals = zeros(num_elec,1);
        vals(n) = 1;

        Svals = EvaluateSplineSurface(vals,Smat);
        %lap_vals = vals-mean(Svals(lapeval_surf_idx),1)'
        for m = 1:num_eval
            M_evalpos(m,n) = Svals(lapeval_center_idx(m))-mean(Svals(lapeval_surround_idx(:,m)),1)';
        end
        
    end
end
