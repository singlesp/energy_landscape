function X = spm_sample_vol_kj(V,xyz,interporder)
X = zeros(size(xyz,2),numel(V));
for i = 1:numel(V)
    %X = spm_sample_vol(V,xyz(:,1),xyz(:,2),xyz(:,3),interporder);
    X(:,i) = spm_sample_vol(V(i),xyz(1,:),xyz(2,:),xyz(3,:),interporder);
end
