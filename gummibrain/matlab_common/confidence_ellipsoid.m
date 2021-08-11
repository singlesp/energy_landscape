function ellipseinfo = confidence_ellipsoid(data,varargin)


p = inputParser;
p.addParamValue('alpha',.95);
p.addParamValue('std',[]);

p.parse(varargin{:});
r = p.Results;

if(~isempty(r.std))
    ci_alpha = 1-(normcdf(-r.std,0,1)*2);
else
    ci_alpha = r.alpha;
end

p=size(data,2); %dimensionality
n=size(data,1);
k=finv(ci_alpha,p,n-p)*p*(n-1)/(n-p);


m1 = mean(data,1);

data1m = data - repmat(m1,size(data,1),1);

[U1,S1,V1] = svd(data1m,0);

%ci_scale = -norminv((1-ci_alpha)/2,0,1);
%pc1 = U1*S1;
%se1 = std(pc1)*ci_scale;

ellipseinfo = [];

for j = 1:numel(k)
    se1=diag(S1)*sqrt(k(j)/(n-1));
    ax1 = diag(se1)*V1.';
    if(p == 2)
        th = linspace(0,2*pi);
        xy = [cos(th(:)) sin(th(:))];
        %el1 = xy*diag(se1)*V1';
        el1 = xy*ax1;
    else
        el1 = [];
    end
    %mean(inpolygon(data1m(:,1),data1m(:,2),el1(:,1),el1(:,2)))

    ei = struct('mean',m1,'alpha',ci_alpha(j),'axpc',ax1,'ellipsexy',el1);
    if(isempty(ellipseinfo))
        ellipseinfo = ei;
    else
        ellipseinfo(end+1) = ei;
    end
end
