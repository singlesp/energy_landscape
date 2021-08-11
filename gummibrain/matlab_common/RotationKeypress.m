%%% switches axis into rotation mode when you press "r" and out of rotation
%%% mode when right click

function RotationKeypress(fig)

if(nargin < 1)
    fig = gcf;
end

try
    id1 = iptaddcallback(fig,'KeyPressFcn',@keypress);

    hrot = rotate3d(fig);
    hCM_rotate3d = uicontextmenu('Callback',@context_callback);
    set(hrot,'uicontextmenu',hCM_rotate3d);

catch
end

%%%% process right clicks while rotating
function context_callback(gcbo,eventdata)
rotate3d off;

%%%% process keypress when not rotating
function keypress(gcbo,eventdata)
if(eventdata.Character == 'r')
    rotate3d on;
end