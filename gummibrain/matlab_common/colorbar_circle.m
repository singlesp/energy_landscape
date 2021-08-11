function hbar = colorbar_circle(cmap,varargin)
p = inputParser;
p.addParamValue('update',[]);
p.addParamValue('parent',gcf);
p.addParamValue('mode','polar');
p.addParamValue('size',101);
p.addParamValue('aa',2);
p.addParamValue('clim',[-pi pi]);
p.addParamValue('start',0);
p.addParamValue('dir',1);
p.addParamValue('ztop',true);
p.addParamValue('colorbar',false);


p.parse(varargin{:});
r = p.Results;

existing_cbar = r.update;
if(~isempty(existing_cbar) && ishandle(existing_cbar))
    r = mergestruct(get(existing_cbar,'userdata'),r);
else
    r.update = [];
    existing_cbar = [];
end
parent = r.parent;
sz = r.size;
aascale = r.aa;
cmap_mode = r.mode;
c_lim = r.clim;
direction = r.dir;
start = r.start;
show_colorbar = r.colorbar;
do_ztop = r.ztop;

if(numel(sz) < 2)
    sz = [sz sz];
end

if(ischar(show_colorbar))
    if(strcmpi(show_colorbar,'on'))
        show_colorbar = true;
    else
        show_colorbar = false;
    end
end

cmapsize = size(cmap,1);
%cmap = circshift(cmap,round([cmapsize/2 0]));
bgcolor = uint8(255*get(gcf,'color'));

szr = 2*floor(sz*aascale/2) + 1;
[x y] = meshgrid(linspace(-1,1,szr(1)),linspace(-1,1,szr(2)));
x = x-mean(x([1 end]));
y = y-mean(y([1 end]));

rho = sqrt(x.^2 + y.^2);
theta = atan2(y,x);


bgimg = repmat(reshape(bgcolor,1,1,3),szr);
switch lower(cmap_mode)
    case {'polar','pol'}
        val = direction*theta + start;
        val = mod(val+pi,2*pi)-pi;
        img = uint8(255*val2rgb(val,cmap,[-pi pi]));
    case {'eccen','eccentricity','ecc'}
        %val = direction*rho/max(rho(:)) + start;
        val = 2*pi*direction*rho + start;
        val = mod(val+pi,2*pi)-pi;
        img = uint8(255*val2rgb(val,cmap,[-pi pi]));
    otherwise
        img = bgimg;
end

circidx = find(repmat(rho,[1 1 3])>1);
img(circidx) = bgimg(circidx);

alphaimg = ones(size(img,1),size(img,2));
alphaimg(rho >= 1) = 0;

img_orig = imresize(img,1/aascale);
alpha_orig = imresize(alphaimg,1/aascale);

fig = parent;
%orig_fig = gcf;
%orig_ax = gca;
%figure(fig);

if(~isempty(existing_cbar))
    hbar = existing_cbar;
    set(fig,'currentaxes',hbar);
    himg = findobj(hbar,'type','image');
    himg = himg(1);
    set(himg,'cdata',img_orig);
else

    hbar = axes('position',[0 0 .2 .2]);
    himg = imshow(img_orig);
end

set(himg,'cdatamapping','scaled','alphadata',alpha_orig,'alphadatamapping','scaled');
set(hbar,'clim',c_lim,'tag','colorbar','userdata',r);

set(fig,'colormap',cmap);
if(isempty(existing_cbar) && show_colorbar)
    colorbar('peer',hbar);
end

%figure(orig_fig);
%axes(orig_ax);

if(isempty(existing_cbar))
    if(do_ztop)
        ch = get(fig,'children');
        chnew = [hbar; setdiff(ch,hbar)];
        if(numel(chnew) == numel(ch))
            set(fig,'children',[hbar; setdiff(ch,hbar)]);
        end
    end
    if(ishandle(himg))
        set(himg,'hittest','on','buttondownfcn',{@imgclick,hbar});
    end
end

if(~ishandle(hbar))
    hbar = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imgclick(h,ev,hbar)
hfcn = get(hbar,'buttondownfcn');
if(~isempty(hfcn))
    hfcn(hbar,ev);
end
