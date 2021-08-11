function M_laplace = montage_laplacian_transfer(elec_pos, elec_tri, max_laplace_dist)
% build surface laplacian transfer matrix for this specific electrode arrangement.
%   * Each ROW of M is the laplacian around a specific electrode
%
% ex: [elec_tri cap_edge] = HeadMeshTriangles(elec_pos);
%     M_mont = montage_laplacian_transfer(elec_pos, elec_tri);
%
% M_mont(po4,:)*data(elec,samples) ==> laplacian centered at po4

%%%%%
%%%%% note: commented-out code helps remove possible "neighbor" vertices
%%%%% that are too far to be considered (ie... edge vertices that "connect"
%%%%% to the other side of the head if naive triangulation was used.
%%%%% To avoid this problem, use HeadMeshTriangles first, which will
%%%%% automatically remove/flag these inappropriate connections
%%%%%


%
% if(nargin < 3 || isempty(max_neighbor))
%     max_neighbor = 7;
% end

if(nargin < 3 || isempty(max_laplace_dist))
    max_laplace_dist = inf;
end

D = distance(elec_pos',elec_pos');
D(D<1e-6) = inf;
elec_spacing = median(min(D,[],2));


max_laplace_dist = max_laplace_dist*elec_spacing;
max_neighbor = 7;
%max_neighbor = inf;

num_elec = size(elec_pos,1);

M_laplace = zeros(num_elec,num_elec);

for n = 1:num_elec
    if(~isempty(elec_tri))
        elec_nei = elec_tri(any(elec_tri == n,2),:);
    else
        elec_nei = 1:num_elec;
    end
    elec_nei = elec_nei(:);
    elec_nei = unique(elec_nei(elec_nei ~= n));

    %only choose neighbors if they are reasonably close
    nei_dist = D(n,elec_nei);
    %nei_dist = distance(elec_pos(n,:)',elec_pos(elec_nei,:)');
    [~, idx] = sort(nei_dist,'ascend');

    idx = idx(1:min(numel(elec_nei),max_neighbor));
    elec_nei = elec_nei(idx);
    nei_dist = nei_dist(idx);
    
    elec_nei = elec_nei(nei_dist <= max_laplace_dist);
    nei_dist = nei_dist(nei_dist <= max_laplace_dist);
     
    %d = 1./(nei_dist/elec_spacing).^2;
    %d = d/sum(d);
    %get fwhm and weight the electrodes by distance instead
    %but then this defeats purpose of simple hjorth
    
    vals = zeros(1,num_elec);
    
    if(isempty(elec_nei))
        warning(['no electrodes meet adjacency parameters for electrode ' num2str(n)]);
        %vals(n) = nan;
        vals(n) = 0;
    else        
        vals(n) = 1;
        vals(elec_nei) = -1/numel(elec_nei);
    end

    M_laplace(n,:) = vals;
end

M_laplace = M_laplace*4/(elec_spacing^2);