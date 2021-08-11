function eloc = headplot3d_setup(subj, bad_elec, which_head, varargin)
% eloc = headplot3d_setup(subj, bad_elec, which_head, [params])

if(nargin < 2)
    bad_elec = [];
end

if(nargin < 3)
    which_head = 'default';
end
which_head = lower(which_head);

datfile = sprintf('%s/%s_128/%s_new.dat',rivalrydir,subj,subj);

missing_elec = {'veo','heo','ekg','emg'};

%%%%% add bad electrodes to list of electrodes to remove
%%%%%% before trying to spherefit/spline the data
if(~isempty(bad_elec))
    if(size(bad_elec,1) > size(bad_elec,2))
        bad_elec = bad_elec';
    end
    b = cellfun(@num2str,num2cell(bad_elec),'UniformOutput',false);
    missing_elec = {missing_elec{:} , b{:}};
end

eloc = electrode_dat2eeglab(datfile,missing_elec,false);
transform = [0 0 0  0 0 0  1 1 1]; %headplot transform (see headplot.m)

switch which_head
    case 'default'
        splinefile = sprintf('%s/%s_128/eeglab_spline_%d.mat',rivalrydir,subj,numel(eloc));
        headplot('setup',eloc,splinefile,'transform',transform,varargin{:}); 
    case {'skin','subject'}
        splinefile = sprintf('%s/%s_128/%s_spline_%d.mat',rivalrydir,subj,subj,numel(eloc));
        skinfile = sprintf('%s/%s_128/%s_skin.mat',rivalrydir,subj,subj);
        headplot('setup',eloc,splinefile,'meshfile',skinfile,'transform',transform,varargin{:}); 
    case 'sphere'
        splinefile = sprintf('%s/%s_128/sphere_spline_%d.mat',rivalrydir,subj,numel(eloc));
        sph_file = sprintf('%s/sphere_surf_simp.mat',rivalrydir);

        S = load(sph_file);
        surf_pos = S.POS;
        
        loc = [[eloc.X]; [eloc.Y]; [eloc.Z]]';
        [c r] = spherefit(loc);
        
        loc = loc - repmat(c',size(loc,1),1); %recentering
        
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
        Sp.ElectrodeNames = num2str(cellfun(@str2num,{eloc.labels}'));
        Sp.indices = 1:numel(eloc);
        %Sp.transform = [0 -10 0 -0.1000 0 -1.6000 1100 1100 1100];
        Sp.transform = transform;
        Sp.comment = '';
        Sp.headplot_version = 2;
        Sp.newElect = loc;
        
        save(splinefile,'-struct','Sp');
        
        S.POS = S.POS * r;
        S.POS(:,3) = S.POS(:,3) + 30;
        skinfile = sprintf('%s/%s_128/%s_sphere.mat',rivalrydir,subj,subj);
        save(skinfile,'-struct','S');
        %headplot('setup',eloc,splinefile,'meshfile',skinfile); 
end
