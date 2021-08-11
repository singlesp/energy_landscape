function c=fastcorr(A,B,dim)

% fastcorr(A,B) 
% Computes Pearson correlation coefficient for
% each pair of corresponding columns in A and B
%
% fastcorr(A,B,DIM)
% Compute corr coef along dimension DIM (default 1 = columns)

if(nargin==2),
   dim = 1;
   if(size(A,1) > 1)
      dim = 1;
   elseif(size(A,2) > 1)
      dim = 2;
   end;
end;

if(~isequal(size(A),size(B)))
    error('Both inputs must be same dimensions');
end

A=demean(A,dim);
B=demean(B,dim);

c=sum(A.*B,dim)./sqrt(sum(A.^2,dim).*sum(B.^2,dim));
