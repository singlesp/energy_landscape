function M = get_realign_parameters(V)
% copied from spm_realign.m in SPM8
%

n = length(V);
M = zeros(n,6);
for j=1:n
    qq     = spm_imatrix(V(j).mat/V(1).mat);
    M(j,:) = qq(1:6);
end
