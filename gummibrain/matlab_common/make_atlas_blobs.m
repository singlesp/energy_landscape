function atlasblobs = make_atlas_blobs(atlasfile,varargin)
% atlasblobs = make_atlas_blobs(atlasfile,'param',value,...)
%
% Required inputs:
% atlasfile = MRI volume file (eg nii.gz) with labeled voxels
%
% Optional inputs:
% backgroundfile = MRI volume file to extract a midsagittal slice from to
%   display in background (must be in same coordinate space as atlasfile)
% roinames = list of ROI names for later reference
% roihemis = list of 'lh' or 'rh' that will be used to decide which views
%   display each ROI on (default = determine automatically from ROI
%   centers)
% volumesmoothing = true/false do we perform a gentle voxel smoothing on 
%   each ROI before extracting surface? (helpful especially for small ROIs,
%   but downside = ROIs get larger/puffier and might overlap)
% surfsmoothing = iterations of surface vertex smoothing after extracting 
%   from volume. Less necessary if volume smoothing performed (0=none)
%
% Output:
% atlasblobs = struct with all the info needed for displaying surface
% rendering later

args = inputParser;
args.addParameter('backgroundfile',[]);
args.addParameter('roinames',[]);
args.addParameter('roihemis',[]);
args.addParameter('volumesmoothing',false);
args.addParameter('surfacesmoothing',0);

args.parse(varargin{:});
args = args.Results;


atlasblobs=struct();

atlasblobs.backgroundslice=[];
atlasblobs.backgroundposition=[];
if(~isempty(args.backgroundfile) && exist(args.backgroundfile))

    [Vbg,~,bgscales,~,~,bghdr] = read_avw_and_header(args.backgroundfile);
    Vbg_x=squeeze(Vbg(round(end/2),:,:));
    
    atlasblobs.backgroundslice=Vbg_x;
    
    bgijk=[size(Vbg,1)/2 1 1;
        size(Vbg,1)/2 size(Vbg,2) 1;
        size(Vbg,1)/2 size(Vbg,2) size(Vbg,3);
        size(Vbg,1)/2 1 size(Vbg,3)];
    bgxyz=affine_transform(bghdr.mat,bgijk);
    
    atlasblobs.backgroundposition=bgxyz;
end


[V,~,scales,~,~,hdr] =read_avw_and_header(atlasfile);

uval=unique(V(V>0));
FV_all={};
fprintf('Generating %d isosurfaces:\n',numel(uval));
for i = 1:numel(uval)
    fprintf('%d ',i);
    if(mod(i,10)==0)
        fprintf('\n');
    end
    Vtmp=V==uval(i);
    %isoval=0.5;
    if(args.volumesmoothing)
        Vtmp=smooth3(Vtmp);
        %isoval=.5*max(Vtmp(:));
    end
    %FV = isosurface(Vtmp,isoval);
    FV = isosurface(Vtmp);
    
    FV.vertices=FV.vertices(:,[2 1 3]); %isosurface seems to swap X and Y for some reason
    
    FV.vertices=affine_transform(hdr.mat,FV.vertices);
    
    if(args.surfacesmoothing>0)
        FV.conn = vertex_neighbours(FV);
        FV.vertices = mesh_smooth_vertices(FV.vertices, FV.conn,[],args.surfacesmoothing);
    end
    FV_all{i}=FV;
end


fprintf('Done!\n');

atlasblobs.FV=FV_all;

atlasblobs.roilabels=uval;
if(~isempty(args.roinames))
    atlasblobs.roinames=args.roinames;
else
    atlasblobs.roinames={};
end

roimean=cellfun(@(x)mean(x.vertices,1),FV_all,'uniformoutput',false);
roimean=cat(1,roimean{:});
atlasblobs.roicenters=roimean;

if(~isempty(args.roihemis))
    atlasblobs.hemisphere=roihemis;
else

    roimin=cellfun(@(x)min(x.vertices,[],1),FV_all,'uniformoutput',false);
    roimin=cat(1,roimin{:});
    roimax=cellfun(@(x)max(x.vertices,[],1),FV_all,'uniformoutput',false);
    roimax=cat(1,roimax{:});
    
    %atlasblobs.hemisphere=repmat({'both'},numel(FV_all),1);
    %atlasblobs.hemisphere(roimax(:,1)<0)={'lh'};
    %atlasblobs.hemisphere(roimin(:,1)>0)={'rh'};
    atlasblobs.hemisphere(roimean(:,1)<0)={'lh'};
    atlasblobs.hemisphere(roimean(:,1)>=0)={'rh'};
    
end