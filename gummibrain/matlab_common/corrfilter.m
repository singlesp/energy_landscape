function cc = corrfilter(A, template)
% normalized template match for columns of 1D signals
% from Yue Wu, Tufts Univ ECE

if(size(A,1) == 1)
    A = A.';
end

template = template(:);
cc = zeros(size(A));
for i = 1:size(A,2)
    cc(:,i) = corrfilter2(A(:,i), template);
end
