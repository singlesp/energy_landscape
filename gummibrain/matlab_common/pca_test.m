close all;
clc;

N1 = 100; %data set sizes
N2 = 100;
d1 = [0 0]; %class positions
d2 = [1 -3]*1;

data1 = [randn(N1,1)+d1(1) randn(N1,1)+d1(2)]; %class1
data2 = [randn(N2,1)+d2(1) randn(N2,1)+d2(2)]; %class2

data1i = 1:N1;
data2i = N1+(1:N2);

data = [data1; data2];

%%%%% calc PCA basis and sort by its variance
sigma = cov(data);
[pca_v, pca_d] = eig(sigma);
[pca_d pca_idx] = sort(diag(pca_d),'descend');
pca_v = pca_v(:,pca_idx);
%v = new basis set (orthonormal so inv(v) = v')
%d = diag = variance accounted for by each column (axis/basis) of v

%%%%%% project data onto new PCA vectors
pca_data = data*pca_v;


%%%%%% visualize the new basis
basis_scale = 1;
m = mean(data);
pca_basis = [];
for n = 1:numel(pca_d)
    b = basis_scale*pca_v(:,n)*pca_d(n);
    pca_basis = [pca_basis; m - b' m + b'];
end

figure;
subplot(2,1,1);
plot(data(data1i,1),data(data1i,2),'b.',data(data2i,1),data(data2i,2),'r.');
hold on;
for n = 1:numel(pca_d)
    plot(pca_basis(n,[1 3]),pca_basis(n,[2 4]),'--black','LineWidth',2);
end
hold off;
axis image;
title('original data with pca bases plotted');

subplot(2,1,2);
plot(pca_data(data1i,1),pca_data(data1i,2),'b.',pca_data(data2i,1),pca_data(data2i,2),'r.');
axis image;
title('data on principle axes');