function retval = display_atlas_blobs(roivals,atlasblobs,varargin)
% img or fig = display_atlas_blobs(roivals,atlasblobs,'param',value,...)
% 
% Required inputs:
% roivals = Rx1 vector with 1 value per ROI in atlasblobs 
% atlasblobs = struct (or struct array) generated from make_atlas_blobs()
% 
% Optional inputs:
% atlasname = name of atlas to select if atlasblobs is an array
% colormap = colormap name or Cx3 [r g b] colormap matrix
% clim = [min max] value range
% render = true/false. true = render all views off screen and return final
%    composed image. False = return visible figure for rotation and viewing 
% roimask = Rx1 vector of true/false to exclude some ROIs if needed
% backgroundimage = true/false to show pre-determined background slice
%    (selected during make_atlas_blobs())
% backgroundcolor = [r g b] background color for figure
% surfacesmoothing = iterations of extra surface smoothing if needed
% crop = true/false to automatically crop rendered images
%
% Output:
% HxWx3 RGB image if render=true
% figure handle if render=false
 
args = inputParser;
args.addParameter('atlasname',[]);
args.addParameter('colormap',[]);
args.addParameter('clim',[]);
args.addParameter('render',false);
args.addParameter('roimask',[]);
args.addParameter('backgroundimage',true);
args.addParameter('backgroundcolor',[]);
args.addParameter('surfacesmoothing',0);
args.addParameter('crop',true);

args.parse(varargin{:});
args = args.Results;

if(~isempty(args.atlasname) && numel(atlasblobs)>1)
    atlasidx=find(strcmpi({atlasblobs.atlasname},args.atlasname));
    if(isempty(atlasidx))
        error('Atlas %d not found: %s\n',args.atlasname);
    end
    
    atlasblobs=atlasblobs(atlasidx);
end

if(isempty(args.colormap))
    args.colormap=get(0,'defaultfigurecolormap');
elseif(ischar(args.colormap))
    args.colormap=evalin('caller',sprintf('%s(256)',args.colormap));
end

if(isempty(args.clim))
    args.clim=[nanmin(roivals(:)) nanmax(roivals(:))];
end

if(isempty(args.backgroundcolor))
    if(args.backgroundimage)
        args.backgroundcolor=[0 0 0];
    else
        args.backgroundcolor=[1 1 1];
    end
end


roicolors = val2rgb(roivals,args.colormap,args.clim);

if(args.render)
    fig=figure('visible','off','color',args.backgroundcolor);
else
    fig=figure('color',args.backgroundcolor);
end
fig.InvertHardcopy = 'off';
ax=axes(fig);
axis(ax,'vis3d','equal','off');
set(ax,'NextPlot','add');


%[p,t]=BuildSphere(2);
%FVsph=struct('vertices',p,'faces',t);
%roiradius=5*ones(size(roixyz,1),1);


tmpd=tempname;
mkdir(tmpd);

hemis={'lh','rh'};

hcam=[];
    
img_all={};
for h = 1:numel(hemis)
    hemi=hemis{h};
    
    
    
    viewnames={'lateral','medial'};
    
    if(strcmpi(hemi,'lh'))
        hemimask=ismember(lower(atlasblobs.hemisphere),{'l','lh','left','both',''});
        viewpoints={[-90 0],[90 0]};
        bglayerpos={200,-200};
    elseif(strcmpi(hemi,'rh'))
        hemimask=ismember(lower(atlasblobs.hemisphere),{'r','rh','right','both',''});
        viewpoints={[90 0],[-90 0]};
        bglayerpos={-200,200};

    end
    
    
    cla(ax);
    for i = 1:numel(roivals)
        if(~isempty(args.roimask) && ~args.roimask(i))
            continue;
        end
        if(~hemimask(i) || isnan(roivals(i)))
            continue;
        end
        FV=atlasblobs.FV{i};
        if(args.surfacesmoothing==0)
            verts_new=FV.vertices;
        else
            if(isfield(FV,'conn'))
                conn=FV.conn;
            else
                conn = vertex_neighbours(FV);
            end
            verts_new = mesh_smooth_vertices(FV.vertices, conn,[],args.surfacesmoothing);
        end
        patch(ax,struct('vertices',verts_new,'faces',FV.faces),'linestyle','none','facecolor',roicolors(i,:));
        
        %verts_new=bsxfun(@plus,FVsph.vertices*roiradius(i),roixyz(i,:));
        %patch(struct('vertices',verts_new,'faces',FVsph.faces),...
        %    'linestyle','none','facecolor',roicolors(i,:));
    end
    material(ax,'dull');
    lighting(ax,'phong');
    
    
    if(~args.render)
        camlight(ax,'headlight');
        RotationHeadlight(fig,true);
    end

    hs=[];
    if(args.backgroundimage)
        %hs=surface(ax,zeros(2,2),'facecolor','texturemap','cdata',atlasblobs.backgroundslice,'linestyle','none');
        %set(hs,'xdata',[0 0]);
        %set(hs,'ydata',atlasblobs.backgroundposition([1 2],2)); %[-126 90]);
        %set(hs,'zdata',reshape(atlasblobs.backgroundposition([1 3; 1 3],3),[2 2])); %[-72 108; -72 108]);
        
        %new style to allow interpolating background slice
        bgsz=size(atlasblobs.backgroundslice);
        bgpos=atlasblobs.backgroundposition;
        hs=surface(ax,zeros(bgsz),'facecolor','interp','cdata',atlasblobs.backgroundslice,'linestyle','none');

        %[bgz,bgy]=meshgrid(linspace(bgpos(1,3),bgpos(3,3),bgsz(2)),linspace(bgpos(1,2),bgpos(2,2),bgsz(1)));
        
        set(hs,'xdata',zeros(bgsz));
        set(hs,'ydata',[-126 90]); %bgy); %[-126 90]);
        set(hs,'zdata',[-72 108; -72 108]); %bgz); %[-72 108; -72 108]);
        
        material(hs,'dull');
        colormap(ax,gray(256));
    end
    

    if(~args.render)
        retval=fig;
        return;
    end
    
    for i = 1:numel(viewpoints)
        view(ax,viewpoints{i});
        if(~isempty(hs))
            set(hs,'xdata',ones(size(get(hs,'xdata')))*bglayerpos{i});
        end
        if(~isempty(hcam))
            delete(hcam);
        end
        hcam=camlight(ax,'headlight');
        
        tmpimgfile=sprintf('%s/atlasblob_%s_%s.png',tmpd,hemi,viewnames{i});
        saveas(fig,tmpimgfile);
        img_all{end+1}=imread(tmpimgfile);
        
    end
end
close(fig);
rmdir(tmpd,'s');

if(args.crop)
    croprect={};
    for i = 1:numel(img_all)
        [~,croprect{i}] = CropBGColor(img_all{i},img_all{1}(1,1,:));
    end
    croprect=cat(1,croprect{:});
    croprect=[min(croprect(:,1:2)) max(croprect(:,3:4))];
    for i = 1:numel(img_all)
        img_all{i}=img_all{i}(croprect(1):croprect(3),croprect(2):croprect(4),:);
    end
end

img_new=[[img_all{1}; img_all{2}] [img_all{3}; img_all{4}]];


retval=img_new;
%imwrite(img_new,sprintf('~/Downloads/%s_allviews.png',whichatlas));