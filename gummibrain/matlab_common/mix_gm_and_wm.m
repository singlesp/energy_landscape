clear all;
close all;
clc;

keepratio = -0.1; % keep 10% of points
triplot = @(tri,loc) trimesh(tri, loc(:,1), loc(:,2), loc(:,3));

bvfile_gm = [researchdir '/v1learning2/Anatomy/sbj1/GM_test/sbj1_tal_LH_GM.srf'];
bvfile_wm = [researchdir '/v1learning2/Anatomy/sbj1/GM_test/sbj1_tal_LH_WM.srf'];

bv_gm = BVQXfile(bvfile_gm);
%gm_tal = bv2tal(bv.VertexCoordinate)';
%gm_tri = bv.TriangleVertex;
[gm_tal gm_tri gm_norms] = SimplifyBVSurf(bv_gm, keepratio);

bv_wm = BVQXfile(bvfile_wm);
%wm_tal = bv2tal(bv.VertexCoordinate)';
%wm_tri = bv.TriangleVertex;
[wm_tal wm_tri wm_norms] = SimplifyBVSurf(bv_wm, keepratio);

%%
mid_tal = zeros(size(wm_tal));

gm_curv = SurfaceCurvature(gm_tri,gm_tal);
wm_curv = SurfaceCurvature(wm_tri,wm_tal);

for n = 1:size(wm_tal,1)
    D = distance(wm_tal(n,:)',gm_tal');
    [minval idx] = min(D);
    mid_tal(n,:) = .5*(wm_tal(n,:)+gm_tal(idx,:));
end



%%
%[mid_tal, mid_tri] = meshcheckrepair(mid_tal, wm_tri);
%mid_curv = SurfaceCurvature(mid_tri,mid_tal);

%h = trisurf(gm_tri,gm_tal(:,1),gm_tal(:,2),gm_tal(:,3),'facecolor',[.8 .8 .8],'edgealpha',0,'facealpha',.5);
normsize = 5;
normstep = 20;

subplot(1,2,1);
h = trisurf(gm_tri,gm_tal(:,1),gm_tal(:,2),gm_tal(:,3),gm_curv,'edgealpha',0,'facealpha',1);
hold on;
% plot3([gm_tal(1:normstep:end,1) gm_tal(1:normstep:end,1)+normsize*gm_norms(1:normstep:end,1)]',...
%     [gm_tal(1:normstep:end,2) gm_tal(1:normstep:end,2)+normsize*gm_norms(1:normstep:end,2)]',...
%     [gm_tal(1:normstep:end,3) gm_tal(1:normstep:end,3)+normsize*gm_norms(1:normstep:end,3)]','-b');
% %hold on;

subplot(1,2,2);
h = trisurf(wm_tri,wm_tal(:,1),wm_tal(:,2),wm_tal(:,3),'facecolor',[.7 1 1],'edgealpha',0);
%h = trisurf(mid_tri,mid_tal(:,1),mid_tal(:,2),mid_tal(:,3),mid_curv,'edgealpha',0);
hold on;
% plot3([wm_tal(1:normstep:end,1) wm_tal(1:normstep:end,1)+normsize*wm_norms(1:normstep:end,1)]',...
%     [wm_tal(1:normstep:end,2) wm_tal(1:normstep:end,2)+normsize*wm_norms(1:normstep:end,2)]',...
%     [wm_tal(1:normstep:end,3) wm_tal(1:normstep:end,3)+normsize*wm_norms(1:normstep:end,3)]','-b');


for s=1:2
    subplot(1,2,s);
    axis vis3d equal;

    material dull;
    shading interp;
    lighting phong;

    axis vis3d equal;

    light('Position',[0 0 1],'Style','infinite');
    light('Position',[0 0 -1],'Style','infinite');
    light('Position',[0 1 0],'Style','infinite');
end

%%

[gm_full gm_full_tri gm_full_norms] = SimplifyBVSurf(bv_gm, -1);
[wm_full wm_full_tri wm_full_norms] = SimplifyBVSurf(bv_wm, -1);

mid_full = zeros(size(wm_tal));

bsz = 1000;
for n = 1:bsz:size(wm_tal,1)
    if(~mod(n-1,100))
        fprintf('%f\n',100*n/size(wm_tal,1));
    end
    
    D = distance(wm_tal(n:n+bsz-1,:)',gm_full');
    [minval minidx] = min(D,[],2);
    
    mid_full(n:n+bsz-1,:) = .5*(wm_tal(n:n+bsz-1,:)+gm_full(minidx,:));
end
%%
FV.vertices = gm_full;
FV.faces = gm_full_tri;

FV_gm=refinepatch(FV);

FV.vertices = wm_simp;
FV.faces = wm_simp_tri;

FV_wm=refinepatch(FV);

%%
FVo.vertices = gm_simp;
%FV.faces=[2 3 4; 4 3 1; 1 2 4; 3 2 1];
FVo.faces = gm_simp_tri;

[FV]=refinepatch(FVo);
%%
minpos = deal(min(gm_tal,[],1));
maxpos = deal(max(gm_tal,[],1));
res = .5;

minpos = floor(minpos/res)*res;
maxpos = ceil(maxpos/res)*res;

%xi = minpos(1):res:maxpos(1);
%yi = minpos(2):res:maxpos(2);
%zi = minpos(3):res:maxpos(3);

[x y z] = meshgrid(minpos(1):res:maxpos(1),...
    minpos(2):res:maxpos(2),...
    minpos(3):res:maxpos(3));

khull_gm = convhulln(gm_tal);
khull_wm = convhulln(wm_tal);

%img_gm=surf2vol(gm_tal,gm_tri,xi,yi,zi);
in_gm = inhull([x(:) y(:) z(:)],gm_tal,khull_gm);
in_wm = inhull([x(:) y(:) z(:)],wm_tal,khull_wm);
%%
% minpos = deal(min(wm_tal,[],1));
% maxpos = deal(max(wm_tal,[],1));
% res = .5;
% 
% minpos = floor(minpos/res)*res;
% maxpos = ceil(maxpos/res)*res;
% 
% xi = minpos(1):res:maxpos(1);
% yi = minpos(2):res:maxpos(2);
% zi = minpos(3):res:maxpos(3);

img_wm=surf2vol(wm_tal,wm_tri,xi,yi,zi);
%khull_gm = convhulln(gm_tal);
%khull_wm = convhulln(wm_tal);

%[x y z] = meshgrid(
%find scalp triangles that fall within this hull
%tri_incap = inhull(sph_loc(sph_tri(:),:),elec_sph,khull);
%%
%[FV]=refinepatch(FV);
FV.vertices = wm_full;
FV.faces = wm_full_tri;

FVm.vertices = mid_full;
FVm.faces = wm_full_tri;

FVtmp = FV;
FVtmp.vertices(:,1) = -FVtmp.vertices(:,1);
%subplot(1,2,1);
patch(FVm,'facecolor',[1 0 0],'edgealpha',0);
hold on;
patch(FVtmp,'facecolor',[1 0 0],'edgealpha',0);

    material dull;
    %shading interp;
    lighting phong;

    axis vis3d equal;

    camlight headlight;
    %light('Position',[0 0 1],'Style','infinite');
    %light('Position',[0 0 -1],'Style','infinite');
    %light('Position',[0 1 0],'Style','infinite');
    RotationHeadlight(gcf,false);
    