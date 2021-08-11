function varargout = AdaptiveRLSfast(L,y,f1,f2,order,fs,varargin)
%
%[Hcn1,Hsn1,Hcn2,Hsn2] = AdaptiveRLSfast(L,y,f1,f2,order,fs)
%adaptive RLS filter
%adapted from Tang & Norcia 1995
% = 200x faster than original
% = 30x faster than AdaptiveRLS2 (prealloc, precalc)
%
% L = 1200
% y = 1 column for each signal
% f1 = 10;
% f2 = 13;
% fs = 1000;
% order = 2;
if(isempty(varargin))
	%lambda = L/(L+1)
    lambda = 1;
    do_snr = false;
    mem = [];
else
    p = inputParser;
    p.addParamValue('lambda',1);
    p.addParamValue('snr',false);
    p.addParamValue('memory',[]);

    p.parse(varargin{:});
    r = p.Results;

    lambda = r.lambda;
    mem = r.memory;
    do_snr = r.snr;
end

if(min(size(y))==1)
    y = y(:);
end

do_calcrms = do_snr || nargout == 5;

if(isempty(L) || isinf(L))
    L = size(y,1);
end

if(~isempty(mem) && mem > 0)
    lambda = mem/(mem+1);
end

if(lambda == 1)
    lambdafilt = ones(L,1)/(L/2);
else
    %lambdafilt = (lambda.^(L:-1:1))/L*2; %oops
    lambdafilt = lambda.^(0:L-1);
end
t = repmat((0:size(y,1)-1)'/fs, 1, size(y,2));


Hcn1 = [];
Hsn1 = [];
Hcn2 = [];
Hsn2 = [];
y_rms = [];


if(~isempty(f1) && f1 > 0)
    %differs from Tang paper.  now h_c = y*cos, h_s = -y*sin
    s = -sin(2*pi*order*f1*t);
    c = cos(2*pi*order*f1*t);
    Hcn1 = filter(lambdafilt,1,y.*c);
    Hsn1 = filter(lambdafilt,1,y.*s);
    if(lambda < 1)
        denom_c = filter(lambdafilt,1,c.^2);
        denom_s = filter(lambdafilt,1,s.^2);
        Hcn1 = Hcn1./denom_c;
        Hsn1 = Hsn1./denom_s;
    end
    Hcn1(isnan(Hcn1)) = 0;
    Hsn1(isnan(Hsn1)) = 0;
end


if(~isempty(f2) && f2 > 0)
    s = -sin(2*pi*order*f2*t);
    c = cos(2*pi*order*f2*t);
    Hcn2 = filter(lambdafilt,1,y.*c);
    Hsn2 = filter(lambdafilt,1,y.*s);
    if(lambda < 1)
        denom_c = filter(lambdafilt,1,c.^2);
        denom_s = filter(lambdafilt,1,s.^2);
        Hcn2 = Hcn2./denom_c;
        Hsn2 = Hsn2./denom_s;
    end
    Hcn2(isnan(Hcn2)) = 0;
    Hsn2(isnan(Hsn2)) = 0;
end

if(do_calcrms)
    y_rms = sqrt(filter(lambdafilt,1,y.^2));
    
    y_med = median(y_rms);
    y_mad = median(abs(y_rms - repmat(y_med,size(y_rms,1),1)));
    y_thresh = y_med - 3*y_mad;
    
%     figure;
%     plot(y_rms);
%     hold on;
%     plot(get(gca,'xlim'),[y_med y_med],':k');
%     plot(get(gca,'xlim'),[y_thresh y_thresh],'r');

    y_rms = max(y_rms,repmat(y_thresh,size(y_rms,1),1));
    y_rms(y_rms < eps) = 1;
    
    
    
	%y_rms = y_rms - repmat(mean(y_rms,1),size(y_rms,1),1);
	%y_rms(y_rms < repmat(std(y_rms,[],1),size(y_rms,1),1)/5) = 0;
    
    
    %y_rms = y_rms - mean(y_rms);
    %y_rms(y_rms < std(y_rms)/5) = 0;
% 
%     figure;
%     subplot(2,1,1);
%     plot(sqrt(Hcn1.^2+Hsn1.^2),'b');
%     hold on;
%     plot(sqrt(Hcn2.^2+Hsn2.^2),'r');
%     plot(y_rms,'k');
%     
%     subplot(2,1,2);
%     
%     plot([henv1],'b');
%     hold on;
%     plot([henv2],'r');
end

if(do_snr)
    if(~isempty(Hcn1))
        Hcn1 = Hcn1./y_rms;
        Hsn1 = Hsn1./y_rms;
    end
    if(~isempty(Hcn2))
        Hcn2 = Hcn2./y_rms;
        Hsn2 = Hsn2./y_rms;
    end
%     
%     subplot(2,1,2);
%     plot(sqrt(Hcn1.^2+Hsn1.^2),'b');
%     hold on;
%     plot(sqrt(Hcn2.^2+Hsn2.^2),'r');

end

if(nargout == 1)
    varargout = {sqrt(Hcn1.^2 + Hsn1.^2)};
elseif(nargout == 2)
    varargout = {Hcn1, Hsn1};
elseif(nargout == 4)
    varargout = {Hcn1, Hsn1, Hcn2, Hsn2};
elseif(nargout == 5)
    varargout = {Hcn1, Hsn1, Hcn2, Hsn2, y_rms};
end

