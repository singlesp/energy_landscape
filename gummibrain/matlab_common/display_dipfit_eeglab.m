function [dippos dipmom dippos_full dipmom_full] = display_dipfit_eeglab(t,EEG)

dipole = EEG.dipfit.model;
model = EEG.dipfit;

dippos = reshape([dipole.posxyz],[3 numel(dipole)])';
dipmom = reshape([dipole.momxyz],[3 numel(dipole)])';
diprv = [dipole.rv]';
dipstrength = sqrt(sum(dipmom.^2,2));
dipstrength = dipstrength.*(1-diprv);

figure;
subplot(4,1,1);
plot3(dippos(:,1),dippos(:,2),dippos(:,3));
axis vis3d equal;
subplot(4,1,2);
plot(t,dipstrength);
title('|d|');
subplot(4,1,3);
plot(t,diprv);
title('rv');
subplot(4,1,4);
plot(t,dipstrength.*(1-diprv));
title('mag adjusted by rv');

set(datacursormode(gcf),'UpdateFcn',@datatip_showindex);

% %if just want to extract peak strength for each time bin
% [dippos_unique idx1 idx2] = unique(dippos,'rows');
% dipstrength_unique = dipstrength(idx1);
% timepoints_unique = t(idx1);

% if want to reflect peak strength and timepoint for each timebin
% problem: what if it is pos1, then pos2, then pos1 again... will merge t=1
% and t=3

timepoints = t;
% [d idx1 idx2] = unique(dippos,'rows');
% for i = 1:numel(idx1)
%     whicht = idx1(i);
%     whichpoints = find(idx2 == whicht);
%     [dipstrength(whichpoints) strengthidx] = max(dipstrength(whichpoints));
%     timepoints(whichpoints) = timepoints(strengthidx);
% %    strengthidx
% end
% timepoints_unique = timepoints(idx1);
% [tnew sortidx] = sort(timepoints_unique);
% idx1 = idx1(sortidx);

dippos_unique = [];
dipstrenth_unique = [];
timepoints_unique = [];
dipmodel_unique = EEG.dipfit.model;

tcount = 0;
curpos = [nan nan nan];
for i = 1:numel(timepoints)
    if(~all(dippos(i,:) == curpos))
        tcount = tcount + 1;
        curpos = dippos(i,:);
        dippos_unique(tcount,:) = curpos;
        
        dipstrength_unique(tcount) = -inf;
        timepoints_unique(tcount) = timepoints(i);
        dipmodel_unique(tcount) = EEG.dipfit.model(i);
    end
    
    dipstr = dipstrength(i,:);
    if(dipstr > dipstrength_unique(tcount))
        dipstrength_unique(tcount) = dipstr;
        timepoints_unique(tcount) = timepoints(i,:);
        dipmodel_unique(tcount) = EEG.dipfit.model(i);
    end
end
dipmodel_unique = dipmodel_unique(1:tcount);

% 
% idx1 = 1:numel(dipstrength);
%dippos_unique = dippos(idx1,:);
%dipstrength_unique = dipstrength(idx1);
%timepoints_unique = timepoints(idx1);

timepoints_unique_str = cellfun(@num2str,num2cell(timepoints_unique),'uniformoutput',false);
%dipsize = 1*dipstrength_unique/max(dipstrength_unique);
dipsize = .5;


EEG2 = EEG;
EEG2.dipfit.model = dipmodel_unique;

% pop_dipplot( EEG2,[] ,'dipolesize',30*dipsize,...
%     'color',1:numel(dippos_unique),...
%     'mri','C:\\Users\\Keith\\MATLAB Toolboxes\\eeglab\\plugins\\dipfit2.2\\standard_BESA\\avg152t1.mat','normlen','on');


%hold on;
%tubeplot(dippos,dipstrength/max(dipstrength),1:numel(dipole),[],'edgealpha',0);
%plotsphere(dippos_unique,dipstrength_unique/max(dipstrength),1:numel(dipstrength_unique),'edgealpha',0);
%tubeplot(dippos_unique,dipstrength_unique/max(dipstrength),1:numel(dipstrength_unique),[],'edgealpha',0);
%axis vis3d equal;

% %%

[sources X Y Z XE YE ZE] = dipplot(EEG2.dipfit.model,'coordformat',EEG.dipfit.coordformat,...
    'dipolesize',30*dipsize,...
    'color',timepoints_unique,...
    'mri','avg152t1.mat','normlen','on');

dipfig = gcf;
l = findobj(gca,'type','line');
ddisk = findobj(l,'marker','.');
dbar = findobj(l,'marker','');
%set(l,'visible','off');
set(ddisk,'visible','off')
hold on;
newdippos = reshape([sources.eleccoord],[3 numel(sources)])';
newdipmom = [XE' YE' ZE'] - [X' Y' Z'];

%plot4(newdippos(:,1),newdippos(:,2),newdippos(:,3),timepoints_unique);
%plot3(newdippos(:,1),newdippos(:,2),newdippos(:,3));
%tubeplot(newdippos,.5,timepoints_unique,[],'edgealpha',0);
tubeplot(newdippos,3*dipstrength_unique/max(dipstrength_unique),timepoints_unique,{'cylinders',false},'edgealpha',0);
%tubeplot(newdippos,dipstrength_unique/max(dipstrength_unique),1:numel(dipstrength_unique),[],'edgealpha',0);
%plotdipole(newdippos,normrows(newdipmom),timepoints_unique,'edgealpha',0)

%text(newdippos(:,1)+3,newdippos(:,2),newdippos(:,3),timepoints_unique_str,'backgroundcolor','w')
%text(newdippos(:,1)+3,newdippos(:,2),newdippos(:,3),timepoints_unique_str,'backgroundcolor','w')
text(newdippos(:,1)+3,newdippos(:,2),newdippos(:,3),timepoints_unique_str,'color','w')
%UniformAxes;
colorbar;

if(nargout > 0)
    dippos = newdippos;
end

if(nargout > 1)
    dipmom = newdipmom;
end

if(nargout > 2)

    [sources_full Xfull Yfull Zfull XEfull YEfull ZEfull] = dipplot(EEG.dipfit.model,'coordformat',EEG.dipfit.coordformat,...
        'dipolesize',30*dipsize,...
        'color',t,...
        'mri','avg152t1.mat','normlen','on');
    %delete(findobj(gcf,'type','line'));
    close(gcf);

    newdippos_full = reshape([sources_full.eleccoord],[3 numel(sources_full)])';
    newdipmom_full = [XEfull' YEfull' ZEfull'] - [Xfull' Yfull' Zfull'];

    dippos_full = newdippos_full;
    
    if(nargout > 3)
        dipmom_full = newdipmom_full;
    end
end
