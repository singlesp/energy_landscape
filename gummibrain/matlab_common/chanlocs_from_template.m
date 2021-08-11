function eloc = chanlocs_from_template(eloc, polartemplate)

eloc_template = readlocs(polartemplate);
template_labels = {eloc_template.labels};

eloc_labels = {eloc.labels};
for i = 1:numel(eloc)
   idx = find(strcmpi(template_labels,eloc_labels{i}));
   if(~isempty(idx))
       eloc(i).radius = eloc_template(idx).radius;
       eloc(i).theta = eloc_template(idx).theta;
       eloc(i).X = eloc_template(idx).X;
       eloc(i).Y = eloc_template(idx).Y;
       eloc(i).Z = eloc_template(idx).Z;
   end
end
