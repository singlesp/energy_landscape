function [H hx] = histN(varargin)
% [H,hx]=histN(X)
% [H,hx]=histN(X,n)
% [H,hx]=histN(X,n1,n2...)
% [H,hx]=histN(X,hx1,hx2,...)
%
% imagesc(hx{1},hx{2},H');

%%%%%%%%%% parse inputs

X = varargin{1};
d = size(X,2);
nx = {};
if(nargin == 1)
    nx = num2cell(20*ones(d,1));
elseif(nargin == 2)
    nx = num2cell(varargin{2}*ones(d,1));
    varargin = varargin(3:end);
elseif(nargin >= d+1)
    nx = varargin(2:d+1);
    varargin = varargin(d+2:end);
else
    error('invalid inputs');
end

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

Xbin = zeros(size(X));
for i = 1:d
    X(X(:,i) < hx{i}(1),i) = hx{i}(1);
    X(X(:,i) > hx{i}(end),i) = hx{i}(end);
    [~, Xbin(:,i)] = histc(X(:,i),hx{i});
end

Xbin = num2cell(Xbin,1);

A = sub2ind(nx,Xbin{:});
a = histc(A,1:size(Hx,1));
H = reshape(a,nx);
