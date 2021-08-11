function varargout = readlocsdat(filename, polartemplate, bad_names, recenter)

if(nargin < 2)
    polartemplate = [];
end

if(nargin < 3 || isempty(bad_names))
    bad_names = {'centroid','ref','ref.'};
end

if(nargin < 4)
    recenter = true;
end

landmark_names = {'nasion','left','right','nz','lpa','rpa'};


M = importdata(filename);


eloc = struct;
%ecount = 0;
for i = 1:size(M.data,1)
    if(i <= numel(M.textdata))
        eloc(i).labels = M.textdata{i};
    else
        eloc(i).labels = '';
        %ecount = ecount + 1;
        %eloc(i).labels = num2str(ecount);
    end
    eloc(i).type = M.data(i,1);
    eloc(i).Y = -M.data(i,2);
    eloc(i).X = M.data(i,3);
    eloc(i).Z = M.data(i,4);
end
eloc = convertlocs(eloc,'cart2all');

[~,landmark_idx] = intersect(upper({eloc.labels}),upper(landmark_names));
trace_idx = find([eloc.type] == 32);
[~,good_idx] = setdiff(upper({eloc.labels}),upper(bad_names));
elec_idx = intersect(find([eloc.type] == 69),good_idx);

landmark_idx = sort(landmark_idx);
trace_idx = sort(trace_idx);
elec_idx = sort(elec_idx);

xyz = [[eloc.X]' [eloc.Y]' [eloc.Z]'];
if(recenter)
    xyz = xyz - repmat(mean(xyz(elec_idx,:),1),size(xyz,1),1);
end

for i = 1:numel(eloc)
    eloc(i).X = xyz(i,1);
    eloc(i).Y = xyz(i,2);
    eloc(i).Z = xyz(i,3);
end

landmarks = eloc(landmark_idx);
eloc = eloc(elec_idx);

if(~isempty(polartemplate))
    eloc_template = readlocs(polartemplate);
    template_labels = {eloc_template.labels};

    eloc_labels = {eloc.labels};
    %for i = 1:numel(eloc)
    %    idx = find(strcmpi(template_labels,eloc(i).labels));
    %    if(~isempty(idx))
    %        eloc(i).radius = eloc_template(idx).radius;
    %        eloc(i).theta = eloc_template(idx).theta;
    %    end
    %end
    
    tempxyz = [[eloc_template.X]' [eloc_template.Y]' [eloc_template.Z]'];
    elecmatch = false(numel(eloc_template),1);
    for i = 1:numel(eloc_template)
        idx = find(strcmpi(eloc_labels,template_labels{i}));
        if(~isempty(idx))
            eloc_template(i).X = eloc(idx).X;
            eloc_template(i).Y = eloc(idx).Y;
            eloc_template(i).Z = eloc(idx).Z;
            elecmatch(i) = true;
        end
    end
    eloc = eloc_template;
    
    if(~all(elecmatch))
        fprintf('\n\n********************************\n');
        fprintf('some electrodes have no digitized position.  fitting template (imprecise):\n');
        fprintf('%s ',template_labels{~elecmatch});
        fprintf('\n********************************\n\n');
        
        tempxyz = [tempxyz ones(size(tempxyz,1),1)]';

        digxyz = [[eloc.X]' [eloc.Y]' [eloc.Z]'];
        digxyz = [digxyz ones(size(digxyz,1),1)]';
        
        badfit = false(size(elecmatch));
        for b = 1:3
            tofit = elecmatch & ~badfit;
            temp2dig_xform = digxyz(:,tofit)/tempxyz(:,tofit);
            tempnew = temp2dig_xform*tempxyz;
            
            d = sum((digxyz(1:3,:)'-tempnew(1:3,:)').^2,2);
            badfit = d > (median(d)+2*1.4826*mad(d,1)); %2s.d. above med
        end
% 
%         figure; 
%         subplot(1,2,1);
%         plot3d(digxyz(1:3,:)','.'); 
%         hold on; 
%         plot3d(tempnew(1:3,:)','rx'); 
%         plot3([digxyz(1,badfit & elecmatch); tempnew(1,badfit & elecmatch)], ...
%             [digxyz(2,badfit & elecmatch); tempnew(2,badfit & elecmatch)], ...
%             [digxyz(3,badfit & elecmatch); tempnew(3,badfit & elecmatch)], 'k');
%         %plot3d(digxyz(1:3,elecmatch)','ko');
%         plot3d(tempnew(1:3,~elecmatch)','ko');
%         axis vis3d equal
%         
%         subplot(1,2,2);
%         [median(d) 2*1.4826*mad(d,1)]
%         hist(d(tofit));
        
        tempnew = tempnew(1:3,:)';
        
        [eloc(~elecmatch).X] = splitvars(tempnew(~elecmatch,1));
        [eloc(~elecmatch).Y] = splitvars(tempnew(~elecmatch,2));
        [eloc(~elecmatch).Z] = splitvars(tempnew(~elecmatch,3));
    end
 
end
    
if(nargout == 1)
    varargout = {eloc};
elseif(nargout == 2)
    varargout = {eloc landmarks};
elseif(nargout == 3)
    varargout = {eloc landmarks xyz(trace_idx,:)};
end
