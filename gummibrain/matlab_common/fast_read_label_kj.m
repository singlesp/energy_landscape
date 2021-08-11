function [labelmask] = fast_read_label_kj(labelfile,numvertex)
%adapted from freesurfer read_label.m

% open it as an ascii file
fid = fopen(labelfile, 'r') ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',labelfile);
  return;
end

fgets(fid) ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',labelfile);
  return;
end

line = fgets(fid) ;
nv = sscanf(line, '%d') ;
l = fscanf(fid, '%d %f %f %f %f\n') ;
l = reshape(l, 5, nv) ;
l = l' ;

fclose(fid) ;

l = l(l(:,1)>=0,[1 end]);
l(:,1) = l(:,1)+1;

if(nargin < 2 || isempty(numvertex))
    numvertex = max(l(:,1));
end

labelmask = zeros(numvertex,1);
labelmask(l(:,1)) = 1;
