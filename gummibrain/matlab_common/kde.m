function [bandwidth,density,xmesh]=kde(data,n,MIN,MAX)
% Reliable and extremely fast kernel density estimator for one-dimensional data;
%        Gaussian kernel is assumed and the bandwidth is chosen automatically;
%        Unlike many other implementations, this one is immune to problems 
%        caused by multimodal densities with widely separated modes (see example). The
%        estimation does not deteriorate for multimodal densities, because we never assume 
%        a parametric model for the data.
% INPUTS:
%     data    - a vector of data from which the density estimate is constructed;
%          n  - the number of mesh points used in the uniform discretization of the
%               interval [MIN, MAX]; n has to be a power of two; if n is not a power of two, then
%               n is rounded up to the next power of two, i.e., n is set to n=2^ceil(log2(n));
%               the default value of n is n=2^12;
%   MIN, MAX  - defines the interval [MIN,MAX] on which the density estimate is constructed;
%               the default values of MIN and MAX are:
%               MIN=min(data)-Range/10 and MAX=max(data)+Range/10, where Range=max(data)-min(data);
% OUTPUTS:
%   bandwidth - the optimal bandwidth (Gaussian kernel assumed);
%     density - column vector of length 'n' with the values of the density
%               estimate at the grid points;
%     xmesh   - the grid over which the density estimate is computed;
%  Reference: Z. I. Botev, J. F. Grotowski and D. P. Kroese
%             "KERNEL DENSITY ESTIMATION VIA DIFFUSION" ,Submitted to the Annals of Statistics, 2009
%    
%
%  Example:
%           data=[randn(100,1);randn(100,1)*2+35 ;randn(100,1)+55];
%           [bandwidth,density,xmesh]=kde(data,2^14,min(data)-5,max(data)+5);
%           plot(xmesh,density)

%  Notes:   If you have a more reliable and accurate one-dimensional kernel density
%           estimation software, please email me at botev@maths.uq.edu.au


data=data(:); %make data a column vector

if nargin<2 % if n is not supplied switch to the default
    n=2^12;
end
n=2^ceil(log2(n)); % round up n to the next power of 2;

if nargin<4 %define the default  interval [MIN,MAX]
    minimum=min(data); maximum=max(data);
    Range=maximum-minimum;
    MIN=minimum-Range/10; MAX=maximum+Range/10;
end
% set up the grid over which the density estimate is computed;
R=MAX-MIN; dx=R/(n-1); xmesh=MIN+[0:dx:R]; N=length(data);
%bin the data uniformly using the grid define above;
initial_data=histc(data,xmesh)/N;
a=dct1d(initial_data); % discrete cosine transform of initial data
% now compute the optimal bandwidth^2 using the referenced method
I=[1:n-1]'.^2; a2=(a(2:end)/2).^2;

% use either fixed_point iteration or fzero to solve the equation
options=optimset('Display','off');
t_star=abs(fzero(@(t)(fixed_point(abs(t))-abs(t)),0,options));
if (isnan(t_star)==1) | (isinf(t_star)==1)
   t_star=fminbnd(@(t)abs(fixed_point(t)-t),0,.1);
end



    function  [out,time,t_p]=fixed_point(t)

        f(5)=2*pi^10*sum(I.^5.*a2.*exp(-I*pi^2*t));
        for s=4:-1:2
            [f(s),t_p]=functional(f(s+1),s);
        end
        time=(2*N*sqrt(pi)*f(2))^(-2/5);
        out=(t-time)/time;

        function [f,t]=functional(df,s)
            K0=prod([1:2:2*s-1])/sqrt(2*pi); 
             const=(1+(1/2)^(s+1/2))/3; 
            t=(2*const*K0/N/df)^(2/(3+2*s));
            f=2*pi^(2*s)*sum(I.^s.*a2.*exp(-I*pi^2*t));
        end

    end

% smooth the discrete cosine transform of initial data using t_star
a_t=a.*exp(-[0:n-1]'.^2*pi^2*t_star/2);
% now apply the inverse discrete cosine transform
if nargout>1
    density=idct1d(a_t)/R;
end

bandwidth=sqrt(t_star)*R; 





end


function data=dct1d(data)
% computes the discrete cosine transform of the column vector data
[nrows,ncols]= size(data);
% Compute weights to multiply DFT coefficients
weight = [1;2*(exp(-i*(1:nrows-1)*pi/(2*nrows))).'];
% Re-order the elements of the columns of x
data = [ data(1:2:end,:); data(end:-2:2,:) ];
% Multiply FFT by weights:
data= real(weight.* fft(data));
end
function out = idct1d(data)
% computes the inverse discrete cosine transform
[nrows,ncols]=size(data);
% Compute weights
weights = nrows*exp(i*(0:nrows-1)*pi/(2*nrows)).';
% Compute x tilde using equation (5.93) in Jain
data = real(ifft(weights.*data));
% Re-order elements of each column according to equations (5.93) and
% (5.94) in Jain
out = zeros(nrows,1);
out(1:2:nrows) = data(1:nrows/2);
out(2:2:nrows) = data(nrows:-1:nrows/2+1);
%   Reference:
%      A. K. Jain, "Fundamentals of Digital Image
%      Processing", pp. 150-153.
end


