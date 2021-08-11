function se = ste(data,dim)

if(nargin < 2)
    dim = 1;
end

n = size(data,1);
se = std(data,0,dim)/sqrt(n);