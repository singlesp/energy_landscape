function CenterText(window, displayinfo, txt, x, y, fontsize, color, fontname)


if(~exist('x','var') || isempty(x))
    x = displayinfo.center(1);
end

if(~exist('y','var') || isempty(y))
    y = displayinfo.center(2);
end

if(~exist('color','var') || isempty(color))
    color = [];
end

if(exist('fontsize','var'))
    Screen(window,'TextSize',fontsize);
end

if(exist('fontname','var'))
    Screen(window,'TextFont',fontname);
end

txtrect = Screen('TextBounds',window,txt);
x = x - txtrect(3)/2;
y = y - txtrect(4)/2;
Screen('DrawText', window, txt, x, y ,color);
