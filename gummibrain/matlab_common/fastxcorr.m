function [c,shifts]=fastxcorr(A,B,dim,maxshift)

% fastcorr(A,B) 
% Computes Pearson correlation coefficient for
% each pair of corresponding columns in A and B
%
% fastcorr(A,B,DIM)
% Compute corr coef along dimension DIM (default 1 = columns)

if(nargin==1),
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

%A=demean(A,dim);
%B=demean(B,dim);

An=sqrt(sum(A.^2,dim));
Bn=sqrt(sum(B.^2,dim));

shifts=-maxshift:maxshift;

c=zeros(numel(shifts),numel(An));

i=0;
for s = shifts
    i=i+1;
    c(i,:)=sum(A.*circshift(B,s,dim),dim)./(An.*Bn);
end
