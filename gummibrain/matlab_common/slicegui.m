function hfig = slicegui(varargin)

if(StringIndex(varargin,'parent'))
    i=StringIndex(varargin,'parent');
    hparent = varargin{i+1};
    varargin = varargin(sort(setdiff(1:numel(varargin),[i i+1])));
else
    hparent = [];
end
gx = [];
gy = [];
gz = [];
if(nargin > 1 && (all(size(varargin{1}) == size(varargin{2})) || min(size(varargin{1})) == 1))
    gx = varargin{1};
    gy = varargin{2};
    gz = varargin{3};
    if(min(size(gx)) == 1)
        if(numel(gx) == 2)
            gx = gx(1):gx(3);
            gy = gy(1):gy(3);
            gz = gz(1):gz(3);
        elseif(numel(gx) == 3)
            gx = gx(1):gx(2):gx(3);
            gy = gy(1):gy(2):gy(3);
            gz = gz(1):gz(2):gz(3);
        end
        [gx gy gz] = meshgrid(gx,gy,gz);
    end
    varargin = varargin(4:end);
end

M = permute(varargin{1},[2 1 3]);
varargin = varargin(2:end);

datasize = size(M);
if(isempty(gx))
    gx_list = 1:datasize(1);
    gy_list = 1:datasize(2);
    gz_list = 1:datasize(3);
else
    xstep = gx(1,2,1)-gx(1,1,1);
    ystep = gy(2,1,1)-gy(1,1,1);
    zstep = gz(1,1,2)-gz(1,1,1);
    gx_list = gx(1,1,1):xstep:gx(1,end,1);
    gy_list = gy(1,1,1):ystep:gy(end,1,1);
    gz_list = gz(1,1,1):zstep:gz(1,1,end);
end

x = gx_list(round(end/2));
y = gy_list(round(end/2));
z = gz_list(round(end/2));
    
havecoords = false;
if(numel(varargin) >= 3 && isnumeric(varargin{1}))
    if(numel(varargin{1}) == 3)
        x = varargin{1}(1);
        y = varargin{1}(2);
        z = varargin{1}(3);
        havecoords = true;
        varargin = varargin(2:end);
    elseif(numel(varargin{1}) == 1)
        x = varargin{1};
        varargin = varargin(2:end);
    end
end
if(~havecoords && numel(varargin) > 0 && isnumeric(varargin{1}))
    y = varargin{1};
    varargin = varargin(2:end);
end
if(~havecoords && numel(varargin) > 0 && isnumeric(varargin{1}))
    z = varargin{1};
    varargin = varargin(2:end);
end

title_str = 'Volume Slice GUI';    
if(isempty(hparent))
    hfig = figure('Name',title_str,'NumberTitle','off','HandleVisibility','on');
else
    hfig = GetParentFigure(hparent);
end

if(isempty(gx))
    hsurf = slice(M,x,y,z);
else
    hsurf = slice(gx,gy,gz,M,x,y,z);
end

for i = 1:numel(hsurf)
    xd = get(hsurf(i),'XData');
    if(all(flatten(xd == xd(1))))
        set(hsurf(i),'userdata',struct('dim',1,'positionfield','XData'));
        continue;
    end
    yd = get(hsurf(i),'YData');
    if(all(flatten(yd == yd(1))))
        set(hsurf(i),'userdata',struct('dim',2,'positionfield','YData'));
        continue;
    end
    set(hsurf(i),'userdata',struct('dim',3,'positionfield','ZData'));
end
axis vis3d equal;
xlabel('x');
ylabel('y');
zlabel('z');

hold on;
hline = [];
hline(1) = plot3(gx_list([1 end]),[y y],[z z],'w');
hline(2) = plot3([x x],gy_list([1 end]),[z z],'w');
hline(3) = plot3([x x],[y y],gz_list([1 end]),'w');

current_point = [x y z];

set(gca,'position',[0 0 1 1]);
current_surface = [];

set(hsurf,'linestyle','none');
if(numel(varargin) > 0)
    set(hsurf,varargin{:});
end
set(hsurf,'ButtonDownFcn',@surf_click);
set(hfig,'WindowButtonUpFcn',@surf_unclick);
set(hfig,'WindowButtonMotionFcn',@fig_drag);
RotationKeypress;

hlabel = uicontrol('parent',hfig,'style','text','units','characters',...
    'position',[0 0 40 3],'string',{'pointer: [0 0 0]','intensity: 0'},'visible','off');
setappdata(hfig,'appdata_slicegui',fillstruct(M,hsurf,current_surface,current_point,...
    datasize,hline,title_str,hlabel,gx_list,gy_list,gz_list));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function surf_click(hObj,event)

hfig = gcbf;
if(strcmpi(get(hfig,'SelectionType'),'alt'))
    return;
end
ax = get(hObj,'parent');

A = getappdata(hfig,'appdata_slicegui');

hsurf = findobj(ax,'type','surface');
set(hfig,'WindowButtonMotionFcn',@fig_drag);
%set(hsurf,'ButtonDownFcn',@surf_click);


p = get(ax,'CurrentPoint');
pfig = get(hfig,'CurrentPoint');
data = get(hObj,'userdata');

surfpos = get(hObj,data.positionfield);
surfpos = surfpos(1);
f = 1-(p(1,data.dim)-surfpos)/(p(1,data.dim)-p(2,data.dim));
px = f*(p(1,:)-p(2,:)) + p(2,:); %point on this specific plane

%px = location in M
%px_mm = location in 3d space
[~,px(1)] = min(abs(A.gx_list - px(1)));
[~,px(2)] = min(abs(A.gy_list - px(2)));
[~,px(3)] = min(abs(A.gz_list - px(3)));
%val =A.M(round(px(2)),round(px(1)),round(px(3)));
val =A.M(round(px(2)),round(px(1)),round(px(3)));

px_mm = [A.gx_list(px(1)) A.gy_list(px(2)) A.gz_list(px(3))];

A.current_surface = hObj;
A.surface_point = px_mm;
A.ax_click_point = p;
A.fig_click_point = pfig;
set(A.hlabel,'string',{sprintf('pointer: [%g %g %g]',px_mm),sprintf('intensity: %.2f',val)},'visible','on');

%t=timer('TimerFcn',{@clearlabel,hfig},'StartDelay',3);
%A.nexttimer = datevec(now);
%A.nexttimer(end) = A.nexttimer(end) + 2.9;
%start(t);
setappdata(hfig,'appdata_slicegui',A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function surf_unclick(hObj,event)
hfig = gcbf;
A = getappdata(hfig,'appdata_slicegui');
A.current_surface = [];
setappdata(hfig,'appdata_slicegui',A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearlabel(hobj,ev,hfig)
A = getappdata(hfig,'appdata_slicegui');
if(etime(datevec(now),A.nexttimer) >= 0)
    set(A.hlabel,'visible','off');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function fig_drag(hObj,event)
A = getappdata(gcbf,'appdata_slicegui');
if(isempty(A.current_surface) || strcmpi(get(gcbf,'SelectionType'),'alt'))
    set(A.hsurf,'ButtonDownFcn',@surf_click);
    return;
end

ax = get(A.current_surface,'parent');
data = get(A.current_surface,'userdata');
surfpos = get(A.current_surface,data.positionfield);

p = get(ax,'CurrentPoint');
pfig= get(hObj,'CurrentPoint');


dfig = pfig-A.fig_click_point;
vec0 = A.ax_click_point(1,:)-A.ax_click_point(2,:);
vec1 = p(1,:)-A.ax_click_point(2,:);
n = zeros(1,3);
n(data.dim) = 1;

c0 = acos((vec0*n')/sqrt(vec0*vec0'));
c1 = acos((vec1*n')/sqrt(vec1*vec1'));

pshift = sign(c0-c1)*sqrt(dfig*dfig');

% px = [];
% [~,px(1)] = min(abs(A.gx_list - A.surface_point(1)));
% [~,px(2)] = min(abs(A.gy_list - A.surface_point(2)));
% [~,px(3)] = min(abs(A.gz_list - A.surface_point(3)));

%newpos = round(A.surface_point(data.dim)+pshift);
newpos_mm = A.surface_point(data.dim)+pshift;

if(data.dim == 1)
    [~,newpos] = min(abs(A.gx_list - newpos_mm));

    newpos = min(max(1,newpos),A.datasize(2));
    newpos_mm = A.gx_list(newpos);
    imgdata = squeeze(A.M(:,newpos,:));
elseif(data.dim == 2)
    [~,newpos] = min(abs(A.gy_list - newpos_mm));
    newpos = min(max(1,newpos),A.datasize(1));
    newpos_mm = A.gy_list(newpos);
    imgdata = squeeze(A.M(newpos,:,:));
elseif(data.dim == 3)
    [~,newpos] = min(abs(A.gz_list - newpos_mm));
    newpos = min(max(1,newpos),A.datasize(3));
    newpos_mm = A.gz_list(newpos);
    imgdata = squeeze(A.M(:,:,newpos));
end

surfpos(:) = newpos_mm;
set(A.current_surface,data.positionfield,surfpos,'CData',imgdata);


A.current_point(data.dim) = newpos_mm;
setappdata(gcbf,'appdata_slicegui',A);

cp = A.current_point;

set(A.hline(1),'ydata',cp([2 2]),'zdata',cp([3 3]));
set(A.hline(2),'xdata',cp([1 1]),'zdata',cp([3 3]));
set(A.hline(3),'xdata',cp([1 1]),'ydata',cp([2 2]));

set(gcbf,'Name',sprintf('%s [%d %d %d]',A.title_str,cp))