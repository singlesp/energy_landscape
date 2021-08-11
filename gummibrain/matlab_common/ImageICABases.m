function bases = ImageICABases(imgs)

%imgs = cell array of images
numimg = numel(imgs);
[imgw imgh] = size(imgs{1});

data = zeros(imgw*imgh,numimg);

for i = 1:numimg
    data(:,i) = reshape(imgs{i},imgw*imgh,1);
end
%for ICA, each COLUMN is a time-point/image/whatever

size(data)

if(0)
%% FastICA
[w s] = runica(data); %check "doc runica" for a description of what w and s are
sep_matrix = w*s;

else
%% RobustICA
tol = 1e-4;     % termination threshold parameter
max_it = 1e3;   % maximum number of iterations per independent component
prewhi = 1;     %prewhitening
dimred = 0;     %dimensionality reduction
deftype = 'r';  %deflation type ('o' orthogonalization, 'r' regression)

%input data = [samples x channels]
[data_ica, H, iter, W] = robustica(data, [], tol, max_it, prewhi, deftype, dimred, [], 1);
%H is the mixing matrix (no seperating matrix given since they already
%separated it for you
sep_matrix = pinv(H);
end

sep_matrix = sep_matrix'; %make it 1 column per component

%normalize each basis (for display mainly)
sep_matrix = sep_matrix - repmat(min(sep_matrix),size(sep_matrix,1),1);
sep_matrix = sep_matrix ./ repmat(max(sep_matrix),size(sep_matrix,1),1);
bases = {};
for i = 1:size(sep_matrix,2)
    bases{i} = reshape(sep_matrix(:,i),imgw,imgh);
end

