function [az_el_tilt xyextent] = generate_sphere_viewpoints(vertsph,roimask)
%   viewroi:        VxR column matrix of R Vx1 binary ROI masks.  If
%                   provided, ignore <viewpts> argument and use
%                   center-of-mass for each ROI as viewpt. (ie: center view
%                   on the ROI)
%   rotateroi:      True|False (default=False). Rotate view so that each
%                   ROI is oriented HORIZONTALLY in the image.  If <tilt>
%                   is given, it is applies AFTER this operation


%need sphere input
%either x,y,z, or ROI, or preset grid info
if(~isempty(options.viewroi))
    viewN=size(options.viewroi,2);
    viewpts=zeros(viewN,3);
    for i = 1:viewN
        viewpts(i,:)=mean(vertsph(options.viewroi(:,i)~=0,:),1);
    end
else
    
if(~isempty(options.viewroi))
    viewN=size(options.viewroi,2);
    viewpts=zeros(viewN,3);
    for i = 1:viewN
        viewpts(i,:)=mean(vertsph(options.viewroi(:,i)~=0,:),1);
    end
else
    
if(~isempty(options.viewroi) && options.rotateroi)
    [roiax,~]=minipca(viewvert(options.viewroi,1:2));
    Mroi=eye(4);
    Mroi(1:2,1:2)=roiax;


    %need to make sure Mroi does not invert vertical axis
    %ie: wrap roitilt to range [-pi/2, pi/2]
    roitilt=asin(roiax(1,2));
    if(acos(roiax(1,1))<0)
        roitilt=pi-roitilt;
    end
    % no mod needed since range is [-pi/2,pi/2]
    % roitilt=mod(roitilt+pi/2,pi)-pi/2;

    Mroi(1:2,1:2)=[cos(roitilt) sin(roitilt); -sin(roitilt) cos(roitilt)];

    viewvert=affine_transform(Mtilt*Mroi,viewvert);
end

if(~isempty(options.viewroi))
    dim2size=max(viewvert(options.viewroi,2))-min(viewvert(options.viewroi,2));
    roipad=1.1;
    dim2size=roipad*dim2size; %padding
    dim1size=sqrt(1.7^2 - dim2size^2);
    newxyextent=[dim1size/roipad dim2size]/2; %not used

end