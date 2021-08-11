function [roc , fdr , auc , auf, tp, fp] = rocanalysis(t, y, n)
%
% ROCANALYSIS : computation of various statistics related to ROC
%
% Copyright 2007, Jean-Philippe Vert
%
% USAGE:
% [roc , fdr , auc , auf] = rocanalysis(t, y, n)
%
% Computes some statistics from a vector or true labals t (+/-1) and of
% predicted values y. n is the size of the vectors to be returned (the roc
% and fdr curves are sampled over 1/n points between 0 and 1)
% roc is the ROC curve (a vector of dimension n), i.e., TP ratio as a
% function of FP ratio. auc is the area under the curve (between 0 and 1)
% fdr is the false discovery rate (between 0 and 1) as a function of the TP
% ratio. It is also a vector of size 1. auf is the area under this curve
% (between 0 and 1).


% process targets

t = t > 0;

% sort by classifier output

[Y,idx] = sort(-y);
t       = t(idx);

% compute true positive and false positive rates

tp = cumsum(t);
fp = cumsum(~t);

% roc

largeroc = zeros(1,fp(end)+1);
largeroc(fp+1)=tp/tp(end);
roc = interp1([0:(1/fp(end)):1] , largeroc , [0:(1/(n-1)):1]);
auc = mean(roc);

% fdr

largefdr = zeros(1,(tp(end)+1));
largefdr(tp(end:-1:1)+1) = fp(end:-1:1) ./ (fp(end:-1:1)+tp(end:-1:1));
fdr = interp1([0:(1/(tp(end))):1] , largefdr , [0:(1/(n-1)):1]);
auf = mean(fdr);

tp = tp/tp(end);
fp = fp/fp(end);