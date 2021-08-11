
function selidx = SelectWeighted(fitvals,n)
%[v,idx] = sort(fitvals,'descend');
if(all(fitvals==fitvals(1)))
    selidx = 1:n;
else
    selidx = randsample(numel(fitvals),n,true,fitvals-min(fitvals));
end