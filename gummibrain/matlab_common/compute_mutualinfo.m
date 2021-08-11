function mi_val = compute_mutualinfo(X_fp,c_fp,discfact)

if(nargin < 3)
    % discretize ecg.data and ecg.icaact as c and X respectively
    discretization_factor = 5;
else
    discretization_factor = discfact;
end
c = convert_int8(round((c_fp-mean(c_fp))*discretization_factor/std(c_fp)));

X = int8(zeros(size(X_fp)));
for i=1:size(X,2)
    X(:,i)=convert_int8(round((X_fp(:,i)-mean(X_fp(:,i)))*discretization_factor/std(X_fp(:,i))));
end
% compute normalized mutual information between each IC and ECG
mi_val = zeros(1,size(X,2));
for i=1:size(X,2)
    mi_val(i)=mutualinfo(c',X(:,i)')/entropy(X(:,i)');
end



function out = convert_int8(in)

in(in>intmax('int8'))=intmax('int8');
in(in<-intmax('int8'))=-intmax('int8');
out=int8(in);
   
