function filenames = freesurfer_write_annot_kj(annotinfo)

filenames = {};
for i = 1:numel(annotinfo)
    a = annotinfo(i);
    verts = (0:numel(a.label_vals)-1)';
    
    if(~all(a.label_rgb(:) == 0) && max(a.label_rgb(:)) <= 1)
        rgb = floor(a.label_rgb*255);
    else
        rgb = floor(a.label_rgb);
    end
    
    label_rgb2int = rgb(:,1) + rgb(:,2)*2^8 + rgb(:,3)*2^16;
    for j = 2:size(rgb,1)
        while(any(label_rgb2int(1:j-1)==label_rgb2int(j)))
            rgb(j,:)=mod(rgb(j,:),255)+1;
            label_rgb2int(j) = rgb(j,1) + rgb(j,2)*2^8 + rgb(j,3)*2^16;
        end
    end
    
    label_rgb2int = rgb(:,1) + rgb(:,2)*2^8 + rgb(:,3)*2^16;
    labels = label_rgb2int(a.label_vals);
    
    ctmat = [rgb zeros(a.numlabels,1) label_rgb2int];
    ct_orig = '';
    ct = struct('numEntries',a.numlabels,...
        'struct_names',{a.label_names},...
        'orig_tab',ct_orig,...
        'table',ctmat);
    
    fname = a.filename;
    write_annotation(fname,verts,labels,ct);
    %atmp = fillstruct(name,numlabels,filename,hemi,label_names,label_vals,label_rgb);
    filenames{end+1} = fname;
end
