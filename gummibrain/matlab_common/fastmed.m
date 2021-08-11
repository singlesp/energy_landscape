function m = fastmed(X,dim)

exact_even = true;

if(nargin < 2)
    dim = 1;
end

if(dim ~= 1)
    sz = size(X);
    sz1 = sz(dim);
    sz = sz([1:dim-1 dim+1:end]);
    X = reshape(X,sz1,[]);
end

m = fastmed_double(double(X),exact_even);
% 
% tic
% nonnans = ~isnan(X);
% nancols = ~all(nonnans,1);
% m = nan(1,size(X,2));
% {'init' toc}
% 
% tic
% m(~nancols) = fastmed_double(X(:,~nancols),exact_even);
% {'nonnans' toc}

%for i = 1:size(X,2)
%    if(nancols(i))
%        m(i) = fastmed_double(double(X(nonnans(:,i),i)),exact_even);
%    end
%end

if(dim ~= 1)
    m = reshape(m,sz);
end