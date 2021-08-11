function [M regressors] = projectionmatrix_kj(datalen,maxpolydeg, covariates)
%this outputs polymatrix*covariatematrix, which is what GLMestimate uses
% to make "combinedmatrix", which is applied to the design matrix
%for the data, GLMdenoise applies covariatematrix*polymatrix*data
%
if(nargin < 2)
    maxpolydeg = 2;
end

if(nargin < 3)
    covariates = [];
end

P = [];
addglmdenoisepath;

%if(maxpolydeg > 0)
    P = constructpolynomialmatrix(datalen,0:maxpolydeg);
%end

if(~isempty(covariates))
    P = [P covariates];
    %C = projectionmatrix(covariates);
    %P = P*C;
end

if(isempty(P))
    M = eye(datalen);
else
    constcol = find(all(P==repmat(P(1,:),size(P,1),1),1));
    if(numel(constcol) > 1)
        P(:,constcol(2:end)) = [];
    end
    %P = spm_orth(P); %this has no effect
    M = projectionmatrix(P);
end

if(nargout > 1)
    regressors = P;
end
rmglmdenoisepath;
