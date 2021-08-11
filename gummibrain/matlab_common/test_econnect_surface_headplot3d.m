clear all;
skinFullFilename = 'italyskin.mat';
%load('/home/kjamison/eConnectome/sample data/EEG/Human/Standard/MAT File/VEP-Standard-Labels.mat');
EEG = pop_matreader('/home/kjamison/eConnectome/sample data/EEG/Human/Standard/MAT File/VEP-Standard-Labels.mat');
load(skinFullFilename);

skindata = struct;
skindata.POS = italyskin.Vertices;
skindata.TRI1 = italyskin.Faces;
output_splinefile = 'test_spline.mat';
output_skinfile = 'test_skin.mat';

model.k = cell2mat({EEG.locations(EEG.vidx).italyskinidx});
elec_xyz = italyskin.Vertices(model.k,:);
	
[c r] = spherefit(elec_xyz);

loc = elec_xyz - repmat(c',size(elec_xyz,1),1); %recentering
%loc = loc * 1.5;
elec_xyz = loc + repmat(c',size(elec_xyz,1),1);

eloc = EEG.locations;
for i = 1:numel(eloc)
	eloc(i).labels = EEG.labels{i};
	eloc(i).X = elec_xyz(i,1);
	eloc(i).Y = elec_xyz(i,2);
	eloc(i).Z = elec_xyz(i,3);
end


		
% newloc = eloc(1);
% newloc

bad_elec = StringIndex({eloc.labels},{'tp9','tp10'});
good_elec = 1:numel(eloc);
good_elec(bad_elec) = [];


[Smat Sinv] = headplot3d_setup(skindata, output_splinefile, output_skinfile, eloc(good_elec));
%%
clc;
figure('color',[1 1 1]);
t = 170;
%t = size(EEG.data,2);
show_labels = false;
do_skinblend = true;
splinevals = [];
existing_ax = [];
cmap = [];
minmax = [];

cmap_levels = 128;
contour_levels = 10;

%minmax = [0 max(EEG.data(good_elec,t))];
minmax = max(max(abs(EEG.data(good_elec,:))));
%cmap = redskinblue(cmap_levels);


[~, Smat] = EvaluateSplineSurface([],output_splinefile);
ax=headplot3d(EEG.data(good_elec,t),output_skinfile,output_splinefile,show_labels,do_skinblend,splinevals,minmax,existing_ax,cmap);%,'SpecularStrength',0.2,'DiffuseStrength',0.8,'AmbientStrength',0);
%hsurf=SetVertexColors(ax,Smat*EEG.data(good_elec,t));
%set(hsurf,'CDataMapping','scaled');

set(findobj(ax,'type','patch'),'SpecularStrength',0.1);%,'DiffuseStrength',0.8);%,'AmbientStrength',0.5);
lighting phong; % phone, gouraud
lightcolor = .25*[0.5 0.5 0.5];
%light('Position',[0 0 1],'color',lightcolor);
light('Position',[0 1 0],'color',lightcolor);
light('Position',[0 -1 0],'color',lightcolor);

surf_vals = Smat*EEG.data(good_elec,t);
levels = linspace(-minmax,minmax,contour_levels);
colors = val2color(levels*1.5,redskinblue(contour_levels));
%colors = [.5 .5 .5];
plot_mesh_contour(skindata.TRI1,skindata.POS,surf_vals,levels,colors,2);
clean_headplot(ax,show_labels,[],[],false);
%clean_headplot(ax, show_electrodes, color_blend_mode, cminmax, reset, cmap)
%%

N = 360;
az = 180+45+linspace(0,360,N+1);
az = az(1:N);
%t = round(linspace(1,size(EEG.data,2),N));
t = linspace(1,size(EEG.data,2),N);
vals = interp1(1:size(EEG.data,2),EEG.data(good_elec,:)',t,'spline')';
splinevals = EvaluateSplineSurface(vals,Smat);

imgfmt = '/home/kjamison/movframe_%04d.png';
%for az = 1:5:360
for i = 1:N
%for i = 250:500

	delete(findobj(ax,'tag','contourline'));
	plot_mesh_contour(skindata.TRI1,skindata.POS,splinevals(:,i),levels,colors,1);
	ax=headplot3d([],[],[],show_labels,do_skinblend,splinevals(:,i),minmax,ax,cmap);
	
	
	%F = TriScatteredInterp(elec_xyz(good_elec,:),vals(:,i),'linear');
	%surf_vals = F(skindata.POS);
	%surf_vals(isnan(surf_vals)) = splinevals(isnan(surf_vals),i);
	%delete(findobj(ax,'tag','contourline'));
	%plot_mesh_contour(skindata.TRI1,skindata.POS,surf_vals,levels,colors,1);
	%ax=headplot3d([],[],[],show_labels,do_skinblend,surf_vals,minmax,ax,cmap);


	%clean_headplot(ax,false,[],[],false);
	view([az(i) 30]);
	
	
	export_fig(sprintf(imgfmt,i),'-a2','-nocrop');
	%img = export_fig(gcf,'-a2','-nocrop');

	%pause(.01);
end

