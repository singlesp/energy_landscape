function eloc = headplot3d_setup(subj, bad_elec, which_head)
% eloc = headplot3d_setup(subj, bad_elec, which_head)

if(nargin < 2)
    bad_elec = [];
end

if(nargin < 3)
    which_head = 'default';
end
which_head = lower(which_head);

datfile = sprintf('%s_128/%s.dat',subj,subj);

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

eloc = electrode_dat2eeglab(datfile,missing_elec);

switch which_head
    case 'default'
        splinefile = sprintf('%s_128/eeglab_spline_%d.mat',subj,numel(eloc));
        headplot('setup',eloc,splinefile); 
    case {'skin','subject'}
        splinefile = sprintf('%s_128/%s_spline_%d.mat',subj,subj,numel(eloc));
        skinfile = sprintf('%s_128/%s_skin.mat',subj,subj);
        headplot('setup',eloc,splinefile,'meshfile',skinfile); 
    case 'sphere'
        splinefile = sprintf('%s_128/sphere_spline_%d.mat',subj,numel(eloc));
        sph_file = 'sphere_surf.mat';

        S = load(sph_file);
        surf_pos = S.POS;
        
        loc = [[eloc.X]; [eloc.Y]; [eloc.Z]]';
        [c r] = spherefit(loc);
        
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
        Sp.transform = [0 -10 0 -0.1000 0 -1.6000 1100 1100 1100];
        Sp.comment = '';
        Sp.headplot_version = 2;
        Sp.newElect = loc;
        
        save(splinefile,'-struct','Sp');
        
        S.POS = S.POS * r;
        S.POS(:,3) = S.POS(:,3) + 30;
        skinfile = sprintf('%s_128/%s_sphere.mat',subj,subj);
        save(skinfile,'-struct','S');
        %headplot('setup',eloc,splinefile,'meshfile',skinfile); 
end
