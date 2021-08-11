% get skin surface
skinFullFilename = 'italyskin.mat';
%load('/home/kjamison/eConnectome/sample data/EEG/Human/Standard/MAT File/VEP-Standard-Labels.mat');
EEG = pop_matreader('/home/kjamison/eConnectome/sample data/EEG/Human/Standard/MAT File/VEP-Standard-Labels.mat');
model.skin = load(skinFullFilename);
%conn =meshconn(model.skin.italyskin.Faces,size(model.skin.italyskin.Vertices,1));

%%
close all;
figure('color',[1 1 1]);
hObject = gcf;
hfig = gcf;
handles.sceneaxes = gca;
options.iselectrode = true;
options.islabel = true;

%EEG = [];


% bad_elec = StringIndex({eloc.labels},{'tp9','tp10'});
% good_elec = 1:numel(eloc);
% good_elec(bad_elec) = [];

% get electrode labels and positions
if isempty(EEG)
    model.electrodes.labels = {};
    model.electrodes.locations = [];
else
    model.k = cell2mat({EEG.locations(EEG.vidx).italyskinidx});
    model.electrodes.labels = EEG.labels(EEG.vidx);
    model.electrodes.locations = model.skin.italyskin.Vertices(model.k,:);
end


axes(handles.sceneaxes);
axcolor = get(hObject, 'color');
set(handles.sceneaxes, 'color', axcolor);
set(handles.sceneaxes, 'DataAspectRatio',[1 1 1]);
model.colorbar = 0;
set(handles.sceneaxes, 'userdata', model);
box off;
axis off;
hold on;
axis vis3d equal;

bad_elec = StringIndex(model.electrodes.labels,{'tp9','tp10'});
good_elec = 1:numel(model.electrodes.labels);
good_elec(bad_elec) = [];

EEG.vidx = good_elec;

if ~isempty(EEG)
    if options.iselectrode
		
        electrcolor = [0.0  0.0  1.0];
		location = 1.1*model.electrodes.locations(EEG.vidx,:);
        plot3(location(:,1), ...
              location(:,2), ... 
              location(:,3), ... 
              'k.','LineWidth',4,'color', electrcolor);
    end

    vnum = length(model.electrodes.labels);

    if options.islabel
        for i = 1:vnum
			if(~any(EEG.vidx == i))
				continue;
			end
            location = 1.05*model.electrodes.locations(i,:);
            text( location(1), location(2), location(3), ... 
                  char(upper(model.electrodes.labels{i})),'FontSize',8 ,...
                  'HorizontalAlignment','center');
        end
    end
else
    vnum = 10;
end
%  
% smoothiter = 5;
% useralpha = .5;
% usermethod = 'laplacian';
% userbeta = .5;
% verts_smooth = smoothsurf(model.skin.italyskin.Vertices,[],conn,smoothiter,useralpha,usermethod,userbeta);
% %verts_smooth = model.skin.italyskin.Vertices;


%%%%%%%%%%%%%%%%
% remove electrodes that are not standard 10-20 ones,
% then get electrode positions, labels and indices on the italyskin.

model.italyskinxy = load('italyskin-in-xy.mat');
model.italyskinxyz = load('italyskin-in-xyz.mat');

model.k = cell2mat({EEG.locations(EEG.vidx).italyskinidx});
model.electrodes.labels = EEG.labels(EEG.vidx);
model.electrodes.locations = model.skin.italyskin.Vertices(model.k,:);

model.X = model.italyskinxy.xy(model.k,1); % standard xy coordinates relative to electrodes on the skin
model.Y = model.italyskinxy.xy(model.k,2);   
zmin = min(model.italyskinxyz.xyz(model.k,3));
zmin = -inf;
Z = model.italyskinxyz.xyz(:,3);
model.interpk = find(Z > zmin); % focus interpolated vertices

model.XI = model.italyskinxy.xy(model.interpk,1);
model.YI = model.italyskinxy.xy(model.interpk,2);
setappdata(hfig,'model',model);

ischanged = 0;
setappdata(hfig,'ischanged',ischanged);

currentpoint = 210;

values = EEG.data(EEG.vidx,currentpoint);
VI = griddata(model.X,model.Y,values,model.XI,model.YI,'v4');

maxV = max(abs(VI(:)));
minV = -maxV;

cmap = colormap;
len = length(cmap);
coef = (len-1)/(maxV - minV);
FaceVertexCData = model.skin.italyskin.FaceVertexCData;
FaceVertexCData(model.interpk,:) = cmap(round(coef*(VI-minV)+1),:);


%%%%%%%%%%
patch('SpecularStrength',0.2,'DiffuseStrength',0.8,'AmbientStrength',0.5,...
     'FaceLighting','phong',...
     'Vertices',model.skin.italyskin.Vertices,...
     'LineStyle','none',...
     'Faces',model.skin.italyskin.Faces,...
     'FaceColor','interp',...
     'FaceAlpha',1,...
     'EdgeColor','none',...
     'FaceVertexCData',FaceVertexCData);
 %%
%clean_headplot(gca,true,'blend',[],true,[])
%%
clc;
figure('color',[1 1 1]);
t = 180;
%t = size(EEG.data,2);
show_labels = false;
do_skinblend = true;
splinevals = [];
existing_ax = [];
cmap = [];
minmax = [];

cmap_levels = 256;
contour_levels = 10;

%minmax = [0 max(EEG.data(good_elec,t))];
minmax = max(max(abs(EEG.data(good_elec,:))));
cmap = redskinblue(cmap_levels);

values = EEG.data(EEG.vidx,t);
%VI = griddata(model.X,model.Y,values,model.XI,model.YI,'v4');

%[~, Smat] = EvaluateSplineSurface([],output_splinefile);
splinevals = VI;
ax=headplot3d([],output_skinfile,output_splinefile,show_labels,do_skinblend,splinevals,minmax,existing_ax,cmap);%,'SpecularStrength',0.2,'DiffuseStrength',0.8,'AmbientStrength',0);


%%

N = 1000;
az = linspace(180,180+360,N);
%t = round(linspace(1,size(EEG.data,2),N));
t = linspace(1,size(EEG.data,2),N);
vals = interp1(1:size(EEG.data,2),EEG.data(good_elec,:)',t,'spline')';
%splinevals = EvaluateSplineSurface(vals,Smat);

splinevals = zeros(size(VI,1),size(vals,2));
for i = 1:size(vals,2)
	splinevals(:,i) = griddata(model.X,model.Y,vals(:,i),model.XI,model.YI,'v4');
end
%%
imgfmt = '/home/kjamison/movframe_%04d.png';
%for az = 1:5:360
%for i = 1:N
for i = 1:270

	delete(findobj(ax,'tag','contourline'));
	plot_mesh_contour(skindata.TRI1,skindata.POS,splinevals(:,i),levels,colors,1);
	ax=headplot3d([],[],[],show_labels,do_skinblend,splinevals(:,i),minmax,ax,cmap);


	%clean_headplot(ax,false,[],[],false);
	%view([az(i) 30]);
	
	
	%export_fig(sprintf(imgfmt,i),'-a2','-nocrop');
	%img = export_fig(gcf,'-a2','-nocrop');

	pause(.01);
end
