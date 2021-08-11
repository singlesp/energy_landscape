% display the brain and allow user to select specific cortex vertices -- 
% click on the brain and the indices are printed out
function SelectLeadfieldVerticesGUI(subject, hemi)

if(nargin < 1)
    subject = 1;
end

if(nargin < 2)
    hemi = 'L';
end

appdata = getappdata(gcf,'appdata');

research_dir = 'C:/Users/Keith/research/';

build_gui = isempty(appdata); %only do this if never drawn the window before
reset_display = false;

%%%% if program just started, or if subject has switched, load data
if(isempty(appdata) || appdata.CurrentSubject ~= subject)
    leadfield_matfile = sprintf(fixslash('%s/currysimulation/sbj%d_pre_leadfield.mat'),research_dir,subject);
    elecpos_file = sprintf(fixslash('%s/currysimulation/diag_horz_avg sbj%d-pre.ele'),research_dir,subject); %exported in TAL
    bv2cdr_file = sprintf(fixslash('%s/v1learning2/Anatomy/sbj%d/sbj%d_BVtoCDR_simp.mat'),research_dir,subject,subject);

    load(bv2cdr_file);
    %bvL_to_cdr Nx1 list of indices for all the leadfield indices in the
    %left hemi of the bv surface
    
    %curryloc 3x20,000; currytri 3x40,000; currylfd 67x60,000
    load(leadfield_matfile);
    curryloc = curryloc';
    currytri = currytri';
    currylfd = currylfd(7:end,:); %skip first 6 rows: 3 landmarks + 3 zeros (?)

    FV = struct;
    FV.vertices = curryloc;
    FV.faces = currytri;
    norms = -patchnormals(FV);

    %%%% leadfield matrix excludes a few points, so remove those
    num_elecs = size(currylfd,1);
    num_points = size(currylfd,2)/3;
    curryloc = curryloc(1:num_points,:);
    currytri = currytri(all(currytri <= num_points,2),:);

    right_hemi_idx = find(curryloc(:,1) < 0);
    right_hemi_idx = bvR_to_cdr;
    right_hemi_tri_idx = false(size(currytri,1),1);
    %right_hemi_tri_idx(bvR_to_cdr) = true;
    for n = 1:numel(right_hemi_idx)
        right_hemi_tri_idx(any(currytri == right_hemi_idx(n),2)) = true;
    end

    %currylfd = reshape(currylfd,[num_elecs 3 num_points]);
    %currylfd_mag = squeeze(sqrt(sum(currylfd.^2,2)));
    
    appdata.CurrentSubject = subject;
    appdata.CurrentHemi = '';
    appdata.right_hemi_tri_idx = right_hemi_tri_idx;
    appdata.curryloc = curryloc;
    appdata.currytri = currytri;
end

if(~strcmpi(appdata.CurrentHemi,hemi))
    if(strcmpi(hemi,'L'))
        goodhemi_tri_idx = ~appdata.right_hemi_tri_idx;
        hemi_sign = -1;
    else
        goodhemi_tri_idx = appdata.right_hemi_tri_idx;
        hemi_sign = 1;
    end
    
    appdata.CurrentHemi = hemi;
    appdata.hemi_sign = hemi_sign;
    appdata.goodhemi_tri_idx = goodhemi_tri_idx;
end


setappdata(gcf,'appdata',appdata);
%%
handles = guidata(gcf);
    
if(build_gui)
    set(gcf,'Renderer','OpenGL');

    %%%%% create black bg for brain image
    img = [0 1; 1 1];
    axbg = axes('units','normalized','position',[0 0 1 1],'visible','off');
    imagesc(img);
    colormap gray;
    uistack(axbg,'bottom');
    axis off;
    axbg2 = axes('units','normalized','position',[0 0 1 1],'visible','off'); %add empty axis to keep bg from rotating
    axis off;

    hAxBrain = axes('units','normalized','position',[0 .5 .5 .5]);
    %hL = trisurf(currytri,curryloc(:,1),curryloc(:,2),curryloc(:,3),currylfd(1,2,:)); 
   
    handles.hAxBrain = hAxBrain;
end

set(gcf,'CurrentAxes',handles.hAxBrain);
hL = trisurf(appdata.currytri(appdata.goodhemi_tri_idx,:),appdata.curryloc(:,1),appdata.curryloc(:,2),appdata.curryloc(:,3));
set(hL,'FaceColor',.8*[1 1 1]);
set(hL,'EdgeAlpha',0);

view([-appdata.hemi_sign*90 0]);  %coronal from rear

%set(gca,'color','k','xtick',[],'ytick',[],'ztick',[]);
camlight headlight;
axis equal vis3d tight;
lighting phong;
material dull;
axis off;

RotationHeadlight(gcf,true);
RotationKeypress(gcf);


global selected_vert_idx;
selected_vert_idx = [];

%selected_vert_idx = [19381 19408 19994 19884]; % for testing purposes

% ch = get(gca,'Children'); 
% chtype = {get(ch(:),'Type')};
% hpatch = ch(strcmpi(chtype{:},'patch'));
set(hL,'ButtonDownFcn',@axclick_surfacenormal);

grid off;

guidata(gcf,handles);
