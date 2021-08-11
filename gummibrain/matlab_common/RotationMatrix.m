function R = RotationMatrix(rx,ry,rz)

cz = cos(rz);
sz = sin(rz);
cx = cos(rx);
sx = sin(rx);
cy = cos(ry);
sy = sin(ry);

Rx = [1 0 0 0;
      0 cx -sx 0;
      0 sx cx 0;
      0 0 0 1];

Ry = [cy 0 -sy 0;
      0 1 0 0;
      sy 0 cy 0;
      0 0 0 1];
  
Rz = [cz -sz 0 0;
      sz cz 0 0;
      0 0 1 0;
      0 0 0 1];
  
R = Rx*Ry*Rz;