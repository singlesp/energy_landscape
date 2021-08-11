function mi = matrix_mutualinfo(X)

redundancy = @(x,y)(2*mutualinfo(x,y)/(entropy(x)+entropy(y)));

mi = zeros(size(X,2));
for i = 1:size(X,2)
    for j = 1:size(X,2)
        %mi(i,j) = mutualinfo(X(:,i),X(:,j));
        mi(i,j) = redundancy(X(:,i),X(:,j));
    end
end
