function [H hx] = ksdensityN(X,ksizes)
% histN(X)
% histN(X,n)
% histN(X,n1,n2...)
% histN(X,hx1,hx2,...)

%%%%%%%%%% parse inputs
if(nargin == 1)
    ksizes = [];
end

d = size(X,2);
nx = num2cell(100*ones(d,1));

hx = {};
for i = 1:d
    if(numel(nx{i}) == 1)
        hx{i} = linspace(min(X(:,i)),max(X(:,i)),nx{i});
    else
        hx{i} = nx{i};
    end
end

%%%%%%%%%%%%%%%%
nx = cellfun(@numel,hx);

Hx = {};
evstr1 = sprintf('Hx{%d},',1:d);
evstr2 = sprintf('1:%d,',nx);
evstr = ['[' evstr1(1:end-1) '] = meshgrid(' evstr2(1:end-1) ');'];
eval(evstr);

Hx = cell2mat(cellfun(@(x)(x(:)),Hx,'uniformoutput',false));

if(isempty(ksizes))
    [f1,xi1,u1] = ksdensity(X(:,1));
    [f2,xi2,u2] = ksdensity(X(:,2));
    
end

Xbin = zeros(size(X));
for i = 1:d
    [~, Xbin(:,i)] = histc(X(:,i),hx{i});
end

Xbin = num2cell(Xbin,1);

A = sub2ind(nx,Xbin{:});
a = histc(A,1:size(Hx,1));
H = reshape(a,nx);
