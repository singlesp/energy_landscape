function [Smat Sinv] = headplot3d_setup(input_skindata, output_splinefile, output_skinfile, eloc, which_surface, varargin)
%[Smat Sinv] = headplot3d_setup(input_skindata, output_splinefile, output_skinfile, eloc, which_surface, varargin)
%
% input_skindata: mesh data or filename.  
%                   a) struct with input_skindata.POS (Nx3) and input_skindata.TRI1 (Tx3) 1-based idx
%               or  b) filename with POS and TRI1 variables defined as
%                       above
%               or blank if using the 'eeglab' surface or 'sphere' 
%
% output_splinefile: file where spline (3d conversion) matrix is saved
%
% output_skinfile: file where mesh file is saved (if you provided data instead of an input filename)
%
% eloc: EEG.chanlocs structure from EEGLAB.  Only include the electrodes you want to plot
%               NOTE: must include 3D X,Y,Z positions!
%
% which_surface: which type of surface do you want to create: 
%               'custom' (default): use subject specific head surface (must include input_skindata)
%               'eeglab': use built-in EEGLAB head (messy)
%               'sphere': interpolate onto a sphere
%
% varargin: any other arguments you might pass to EEGLAB headplot in setup mode
%
% %%%%%%%%%%%%%%%%%%%%%%
% Returns:
%   Smat: topo->spline transform (verts x elecs)
%   Sinv (optional): spline->topo transform (elecs x verts)

if(ischar(input_skindata))
    input_skinfile = input_skindata;
elseif(~isempty(input_skindata))
    M.POS = input_skindata.POS;
    M.TRI1 = input_skindata.TRI1;
    if(min(M.TRI1(:)) == 0)
        M.TRI1 = M.TRI1 + 1;
    end
    save(output_skinfile,'-struct','M');
    input_skinfile = output_skinfile;
else
    input_skinfile = [];
end

if(nargin < 5 || isempty(which_surface))
    if(isempty(input_skinfile))
        which_surface = 'eeglab';
    else
        which_surface = 'custom';
    end
end
which_surface = lower(which_surface);

transform = [0 0 0  0 0 0  1 1 1]; %headplot transform (see headplot.m)


switch which_surface
    case 'eeglab'
        loc = [[eloc.X]; [eloc.Y]; [eloc.Z]]';
        [c r] = spherefit(loc);

        loc = loc - repmat(c',size(loc,1),1); %recentering
        locnorm = loc ./ repmat(sqrt(sum(loc.^2,2)),1,3);
        loc = locnorm * 20;
        loc(:,3) = loc(:,3) + 30;
        
        for i = 1:numel(eloc)
            eloc(i).X = loc(i,1);
            eloc(i).Y = loc(i,2);
            eloc(i).Z = loc(i,3);
        end
        
        headplot('setup',eloc,output_splinefile,'transform',transform,varargin{:}); 
    case {'custom'}
        headplot('setup',eloc,output_splinefile,'meshfile',input_skinfile,'transform',transform,varargin{:}); 
    case 'sphere'
        %sph_file = sprintf('%s/sphere_surf_simp.mat',rivalrydir);
        sph_file = 'sphere_surf_simp.mat';

        S = load(sph_file);
        surf_pos = S.POS;
        
        loc = [[eloc.X]; [eloc.Y]; [eloc.Z]]';
        [c r] = spherefit(loc);
        
        loc = loc - repmat(c',size(loc,1),1); %recentering
        loc = loc * 1.1;
        
        %S.POS = S.POS*r;
        %S.POS(:,3) = S.POS(:,3);
        S.center = [0 0 0];
       
        locnorm = loc ./ repmat(sqrt(sum(loc.^2,2)),1,3);
        
        Sp = struct;
        [Sp.G Sp.gx] = spherespline(surf_pos(:,1),surf_pos(:,2),surf_pos(:,3),...
            locnorm(:,1),locnorm(:,2),locnorm(:,3));
        
        loc(:,3) = loc(:,3) + 30;
        Sp.Xe = loc(:,1);
        Sp.Ye = loc(:,2);
        Sp.Ze = loc(:,3);
        if(isnumeric(eloc(1).labels))
            Sp.ElectrodeNames = num2str(cellfun(@str2num,{eloc.labels}'));
        else
            Sp.ElectrodeNames = cell2mat(cellfun(@(x)(sprintf('% 4s',x)),{eloc.labels}','uniformoutput',false));
        end
        Sp.indices = 1:numel(eloc);
        %Sp.transform = [0 -10 0 -0.1000 0 -1.6000 1100 1100 1100];
        Sp.transform = transform;
        Sp.comment = '';
        Sp.headplot_version = 2;
        Sp.newElect = loc;
        
        save(output_splinefile,'-struct','Sp');
        
        S.POS = S.POS * r;
        S.POS(:,3) = S.POS(:,3) + 30;
        save(output_skinfile,'-struct','S');
        %headplot('setup',eloc,splinefile,'meshfile',skinfile); 
end

if(nargout == 1)
    [~, Smat, ~] = EvaluateSplineSurface([],output_splinefile);
elseif(nargout == 2)
    [~, Smat, Sinv] = EvaluateSplineSurface([],output_splinefile);
end