function bases = ImagePCABases(imgs)

%imgs = cell array of images
numimg = numel(imgs);
[imgw imgh] = size(imgs{1});

data = zeros(numimg,imgw*imgh);

for i = 1:numimg
    data(i,:) = reshape(imgs{i},1,imgw*imgh);
end
   
size(data)

%%%%% calc PCA basis and sort by its variance
sigma = cov(data);
[pca_v, pca_d] = eig(sigma);
[pca_d pca_idx] = sort(diag(pca_d),'descend');
pca_v = pca_v(:,pca_idx);
%v = new basis set (orthonormal so inv(v) = v')
%d = diag = variance accounted for by each column (axis/basis) of v

%normalize each basis (for display mainly)
pca_v = pca_v - repmat(min(pca_v),size(pca_v,1),1);
pca_v = pca_v ./ repmat(max(pca_v),size(pca_v,1),1);
bases = {};
for i = 1:size(pca_v,2)
    bases{i} = reshape(pca_v(:,i),imgw,imgh);
end