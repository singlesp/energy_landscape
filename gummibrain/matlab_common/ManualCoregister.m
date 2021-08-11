function surflist = ManualCoregister(surflist, starting_surface)
%   surflist = ManualCoregister(surflist, starting_surface)
%
%   Object selection:
%       comma key (,): modify previous object
%       period key (.): modify next object
%
%   Translation:
%       Left/Right arrows: move electrodes left/right (x-axis)
%       Up/Down arrows: move electrodes forward/backward (y-axis)
%       Ctrl+Up/Down: move electrodes up/down (z-axis)
%       Shift+[keys]: fine movements
%
%   Rotation:
%       Alt+Left/Right arrows: rotate around y-axis (tilt left/right)
%       Alt+Up/Down arrows: rotate around x-axis (tilt forward/backward)
%       Alt+Ctrl+Up/Down: rotate around z-axis
%       Shift+[keys]: fine rotation
%
%   Scaling:
%       plus key (+): increase size
%       minus key (-): decrease size
%       zero key (0): reset scale
%       Shift+[keys]: fine scaling

if(nargin < 2)
    current_surface = 1;
else
    current_surface = starting_surface;
end

eval(['help ' mfilename]);


fig = figure('color',[1 1 1]);
%set(fig,'Renderer','OpenGL');
%opengl hardware;
%set(fig,'Renderer','zbuffer')
set(fig,'CloseRequestFcn',@closefcn)
set(fig,'KeyPressFcn',@keypress);

surf_params = {'linestyle','none','facealpha',.5,'edgecolor','k'};
point_params = {'marker','o','color','k'};
text_params = {'color','k'};

surf_params_selected = {'linestyle','-','edgecolor','r'};
point_params_selected = {'color','r'};
text_params_selected = {'color','r'};

default_transform = [1    0    0     0
    0    1   0     0
    0    0    1    0
    0    0    0    1];

try

    hsurf = {};
    for i = 1:numel(surflist)
        if(~isfield(surflist,'transform') || isempty(surflist(i).transform))
            surflist(i).transform = default_transform;
        end
        surflist(i).newverts = affine_transform(surflist(i).transform,surflist(i).verts);
        
        if(isempty(surflist(i).params))
            surflist(i).params = {};
        end
        if(strcmpi(surflist(i).type,'surface'))
            hsurf{i} = trisurf(surflist(i).tri,surflist(i).verts(:,1),surflist(i).verts(:,2),surflist(i).verts(:,3),zeros(size(surflist(i).verts,1),1),'tag',sprintf('surf%d',i),surf_params{:},surflist(i).params{:});
        elseif(strcmpi(surflist(i).type, 'points'))
            hsurf{i} = plot3(surflist(i).verts(:,1),surflist(i).verts(:,2),surflist(i).verts(:,3),'o','tag',sprintf('surf%d',i),point_params{:},surflist(i).params{:});
            if(~isempty(surflist(i).labels))
                %insert labels backwards because that's the order they will appear in the figure
                %children
                text(surflist(i).verts(end:-1:1,1),surflist(i).verts(end:-1:1,2),surflist(i).verts(end:-1:1,3),surflist(i).labels(end:-1:1),'tag',sprintf('surf%d-text',i),text_params{:});
            end
        end

        hold on;
    end
    axis vis3d equal;

    hold off;

    xlabel('x');
    ylabel('y');
    zlabel('z');

    %hold on;
    material dull;
    lighting phong;
    axis equal vis3d;



    %camlight headlight;
    light('Position',[0 0 1],'Style','infinite');
    light('Position',[0 -1 0],'Style','infinite');
    light('Position',[0 1 0],'Style','infinite');

    fprintf('Current object: %d of %d\n',current_surface,numel(surflist));

    A = struct();
    A.axis = gca;
    A.surflist = surflist;
    A.surfchanged = true(size(surflist));
    A.surf_params_default = surf_params;
    A.point_params_default = point_params;
    A.text_params_default = text_params;
    A.surf_params_selected = surf_params_selected;
    A.point_params_selected = point_params_selected;
    A.text_params_selected = text_params_selected;
    A.complete = false;
    A.current_surface = current_surface;
    A.big_shift = 10;
    A.small_shift = 1;
    A.big_rotation = 10 * (pi/180); %10 degrees
    A.small_rotation = .25 * (pi/180); %1 degree
    A.big_scale = 1.1;
    A.small_scale = 1.01;

    setappdata(fig,'appdata',A);

    %RotationKeypress;
    refresh_gui(fig);
    
    waitfor(fig,'userdata');
catch
	surflist = [];
	delete(fig);
	rethrow(lasterror);
	return;

end

A = getappdata(fig,'appdata');
surflist = A.surflist;
delete(fig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keypress(src,event)

    %ignore modifier keypresses (only process them as modifiers)
    if(~isempty(event.Modifier) && any(strcmpi(event.Key,event.Modifier)))
        return;
    end

    A = getappdata(src,'appdata');

    curidx = A.current_surface;
    newidx = curidx;
    
    tra =  A.surflist(curidx).transform(1:3,4);
    rot = A.surflist(curidx).transform(1:3,1:3);
    scale = 1;
    
    is_shift = any(strcmpi('shift',event.Modifier));
    is_ctrl = any(strcmpi('control',event.Modifier));
    is_alt = any(strcmpi('alt',event.Modifier));

    if(is_shift)
        shift_amt = A.small_shift;
        rot_amt = A.small_rotation;
        scale_amt = A.small_scale;
    else
        shift_amt = A.big_shift;
        rot_amt = A.big_rotation;
        scale_amt = A.big_scale;
    end

    c = cos(rot_amt);
    s = sin(rot_amt);

    rot_pos = [c -s; s c]; %positive rotation
    rot_neg = [c s; -s c]; %negative rotation;

    rot_px = eye(3); rot_nx = rot_px; rot_py = rot_px; rot_ny = rot_px; rot_pz = rot_px; rot_nz = rot_px;

    rot_py(2:3,2:3) = rot_pos;
    rot_ny(2:3,2:3) = rot_neg;

    rot_px([1 3],[1 3]) = rot_pos;
    rot_nx([1 3],[1 3]) = rot_neg;

    rot_pz(1:2,1:2) = rot_pos;
    rot_nz(1:2,1:2) = rot_neg;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% scaling %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%
    switch event.Key
        case 'equal'
            scale = scale * scale_amt;
        case 'hyphen'
            scale = scale / scale_amt;
        case '0'
            scale = 1;
        case 'leftbracket'
            %skinalpha = A.skinalpha - .05;
        case 'rightbracket'
            %skinalpha = A.skinalpha + .05;
        case 'period'
            newidx = curidx + 1;
        case 'comma'
            newidx = curidx - 1;
    end
    
    if(newidx ~= curidx)
        newidx = mod(newidx-1,numel(A.surflist))+1;
        A.surfchanged([newidx A.current_surface]) = true;
        A.current_surface = newidx;
        fprintf('Current object: %d of %d\n',A.current_surface,numel(A.surflist));
        setappdata(src,'appdata',A);
        refresh_gui(src);
        return;
    end
    %skinalpha = min(max(skinalpha,0),1);
    %if just a skin alpha change, don't update the rest of the transform
    %if(skinalpha ~= A.skinalpha)
    %    A.skinalpha = skinalpha;
    %    skinsurf = findobj(src,'tag','skinsurf');
    %    set(skinsurf,'FaceAlpha',skinalpha);
    %    setappdata(src,'appdata',A);
    %    return;
    %end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% rotating %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%
    if(is_alt)
        if(is_ctrl) %U/D = +x/-x, L/R = +z/-z
            switch event.Key
                case 'leftarrow'
                    rot = rot_px*rot;
                case 'rightarrow'
                    rot = rot_nx*rot;
                case 'downarrow'
                    rot = rot_pz*rot;
                case 'uparrow'
                    rot = rot_nz*rot;
            end
        else %U/D = +x/-x, L/R = +y/-y
            switch event.Key
                case 'leftarrow'
                    rot = rot_px*rot;
                case 'rightarrow'
                    rot = rot_nx*rot;
                case 'downarrow'
                    rot = rot_py*rot;
                case 'uparrow'
                    rot = rot_ny*rot;
            end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% shifting %%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%
    else 

        if(is_ctrl) %shifting in X-Z plane
            switch event.Key
                case 'leftarrow'
                    tra(1) = tra(1) - shift_amt;
                case 'rightarrow'
                    tra(1) = tra(1) + shift_amt;
                case 'downarrow'
                    tra(3) = tra(3) - shift_amt;
                case 'uparrow'
                    tra(3) = tra(3) + shift_amt;
            end

        else %shifting in X-Y plane
            switch event.Key
                case 'leftarrow'
                    tra(1) = tra(1) - shift_amt;
                case 'rightarrow'
                    tra(1) = tra(1) + shift_amt;
                case 'downarrow'
                    tra(2) = tra(2) - shift_amt;
                case 'uparrow'
                    tra(2) = tra(2) + shift_amt;
            end
        end
    end

    A.surflist(curidx).transform = eye(4);
    A.surflist(curidx).transform(1:3,1:3) = scale*rot;
    A.surflist(curidx).transform(1:3,4) = tra;

    A.surflist(curidx).newverts = affine_transform(A.surflist(curidx).transform, A.surflist(curidx).verts);
    A.surfchanged(curidx) = true;

%     fprintf('Current object: %d of %d\n',A.current_surface,numel(A.surflist));
%     fprintf('\tTranslation: [%.2f %.2f %.2f]\n',tra);
%     fprintf('\tScale: ');
%     fprintf('%.2f ',element(colnorm(A.surflist(curidx).transform(1:3,1:3)),1));
%     fprintf('\n\tRotation: \n');
%     fprintf('\t\t%.2f\t%.2f\t%.2f\n',normcols(A.surflist(curidx).transform(1:3,1:3))');
%     %fprintf('Transform:\n');
%     %fprintf('\t%.2f\t%.2f\t%.2f\t%.2f\n',A.surflist(curidx).transform');
%     fprintf('\n');
%     
  
    %suptitle(src,sprintf('Current object: %d of %d\n',A.current_surface,numel(A.surflist)));
  
    setappdata(src,'appdata',A);
    refresh_gui(src);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function refresh_gui(fig)

    A = getappdata(fig,'appdata');
    surflist = A.surflist;
    
    for i = 1:numel(surflist)
        if(A.surfchanged(i))
            A.surfchanged(i) = false;
        else
            continue;
        end

        hsurf = findobj(fig,'tag',sprintf('surf%d',i));
        switch surflist(i).type
            case 'surface'
                if(i == A.current_surface)
                    params = A.surf_params_selected;
                else
                    params = A.surflist(i).params;
                end
                set(hsurf,'Vertices',surflist(i).newverts,A.surf_params_default{:}, params{:});
            case 'points'
                if(i == A.current_surface)
                    params = A.point_params_selected;
                else
                    params = A.surflist(i).params;
                end
                set(hsurf,'XData',surflist(i).newverts(:,1),'YData',...
                   surflist(i).newverts(:,2),...
                   'ZData',surflist(i).newverts(:,3),A.point_params_default{:},params{:});
        end

        if(~isempty(surflist(i).labels))
            htext = findobj(fig,'tag',sprintf('surf%d-text',i));
            if(i == A.current_surface)
                params = A.text_params_selected;
            else
                params = {};
            end
            for t = 1:numel(htext)
                set(htext(t),'Position', surflist(i).newverts(t,:),A.text_params_default{:}, params{:});
            end
        end
    end
    axis tight;
    setappdata(fig,'appdata',A);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closefcn(src,eventdata)
set(src,'userdata','complete');
