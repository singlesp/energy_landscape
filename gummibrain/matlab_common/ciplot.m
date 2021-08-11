function ciplot(lower,upper,x,color,alpha,edgestyle)
     
% ciplot(lower,upper)       
% ciplot(lower,upper,x)
% ciplot(lower,upper,x,colour)
% ciplot(lower,upper,x,colour,alpha)
% ciplot(lower,upper,x,colour,alpha,edgestyle)
%
% colour can be letter, 3x1 color, or 4x1 color+alpha
% edgestyle can be -, :, --, etc.
%
% Plots a shaded region on a graph between specified lower and upper confidence intervals (L and U).
% l and u must be vectors of the same length.
% Uses the 'fill' function, not 'area'. Therefore multiple shaded plots
% can be overlayed without a problem. Make them transparent for total visibility.
% x data can be specified, otherwise plots against index values.
% colour can be specified (eg 'k'). Defaults to blue.

% Raymond Reynolds 24/11/06
% - modified by Keith Jamison to add style/alpha options

if length(lower)~=length(upper)
    error('lower and upper vectors must be same length')
end

if(nargin<6)
    edgestyle = 'none';
end

if(nargin<5)
    alpha = 1;
end

if nargin<4
    color='b';
elseif numel(color) == 4
    alpha = alpha*color(4);
    color = color(1:3);
end

if nargin<3
    x=1:length(lower);
end

% convert to row vectors so fliplr can work
if find(size(x)==(max(size(x))))<2
x=x'; end
if find(size(lower)==(max(size(lower))))<2
lower=lower'; end
if find(size(upper)==(max(size(upper))))<2
upper=upper'; end

fill([x fliplr(x)],[upper fliplr(lower)],color,'LineStyle',edgestyle,'FaceAlpha',alpha)


