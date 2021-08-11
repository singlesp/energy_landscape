function [p,x]=colorhist(varargin)

symmarg=cellfun(@(x)(ischar(x) && regexpimatch(x,'^symm')),varargin);
is_symmetric=any(symmarg);

varargin=varargin(~symmarg);

[H,hx]=hist(varargin{:});

x=[hx; hx];
x=x(2:end);
y=[H;H];
y=y(1:end-1);
y(end)=0;

x=[x x(end:-1:1)];
y=[y zeros(size(y))];

x=x(:);
y=y(:);

p=patch(x(:),y(:),x(:),'linestyle','-');

c_lim=[min(x) max(x)];
if(is_symmetric)
    c_lim=[-1 1]*max(x);
    % set(gca,'clim',[-1 1]*max(x));
end
set(gca,'clim',c_lim);
