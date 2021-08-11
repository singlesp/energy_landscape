function rms = rms_deviation(affmat1,affmat2,center,rmax)
%implementing rms_deviation from
%/opt/local/fsl-5.0.6/fsl/src/miscmaths/miscmaths.cc
%
isodiff=affmat1*inv(affmat2)-eye(4);
adiff=isodiff(1:3,1:3);
tr=isodiff(1:3,4)+adiff*center(:);
rms = sqrt(tr.'*tr + (rmax*rmax/5)*trace(adiff.'*adiff));