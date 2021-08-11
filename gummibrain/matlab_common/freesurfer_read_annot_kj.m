function annotinfo = freesurfer_parse_annot_kj(annotfiles)

%http://surfer.nmr.mgh.harvard.edu/fswiki/CorticalParcellation
%aparc.annot = Desikan-Killiany Atlas (2006)
% It is a gyral based atlas: ie, a gyrus was defined as running between the bottoms 
%   of two adjacent sulci. That is, a gyrus includes the part visible on the 
%   pial view + adjacent banks of the sulci limiting this gyrus.
%
%aparc.a2009s.annot = Destrieux Atlas (2009)
% It is based on a parcellation scheme that first divided the cortex into 
%   gyral and sulcal regions, the limit between both being given by the curvature 
%   value of the surface. A gyrus only includes the cortex visible on the pial view, 
%   the hidden cortex (banks of sulci) are marked sulcus
%
%aparc.DKTatlas40 = Desikan–Killiany–Tourville (2012)
%
%BA.annot = Brodmann (with V1 = BA17, V2=BA18, MT=BA?)
%
%thresh = ?
%

if(ischar(annotfiles))
    annotfiles = {annotfiles};
end

annotinfo = [];
for i = 1:numel(annotfiles)
    [vidx, label, ctable] = read_annotation(annotfiles{i},0);
    vidx = vidx + 1;
    [u,~,labelidx] = unique(label);
    
    table_idx = zeros(numel(u),1);
    for j = 1:numel(u)
        tidx = find(ctable.table(:,5)==u(j));
        if(~isempty(tidx))
            table_idx(j) = tidx(1);
        end
    end
    
        
    %label_vals = zeros(obj.numverts,1);
    %label_vals(vidx) = labelidx;
    %label_vals = labelidx;
    label_vals = labelidx;
    
    fn = justfilename(annotfiles{i},true);
    [~,~,~,~,m] = regexp(fn,'^([rl]h)\.');
    if(isempty(m))
        hemi = '';
    else
        hemi = m{1}{1};
    end
    name = regexprep(fn,['^' hemi '\.'],'');
    filename = annotfiles{i};
    
    lblnames_temp = ctable.struct_names(table_idx(table_idx>0));
    lblrgb_temp = ctable.table(table_idx(table_idx>0),1:3);
    
    label_names = repmat({'unknown'},numel(table_idx),1);
    label_rgb = repmat([0 0 0],numel(table_idx),1);
    label_names(table_idx>0) = lblnames_temp;
    label_rgb(table_idx>0,:) = lblrgb_temp;
    
    %label_names = ctable.struct_names(table_idx);
    %label_rgb = ctable.table(table_idx,1:3);
    numlabels = numel(table_idx);
    atmp = fillstruct(name,numlabels,filename,hemi,label_names,label_vals,label_rgb);
    annotinfo = [annotinfo atmp];
end
