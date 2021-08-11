function aff = construct_rotmat_euler(params, centre)
%eddy_parameters file is [tx ty tz rx ry rz]
%this function wants: [rx ry rz tx ty tz]
n=numel(params);

aff=eye(4);
if(n <= 0 )
    return;
end
angl=[params(1) 0 0].';
aff=aff*make_rot(angl,centre);

angl=[0 params(2) 0].';
aff=aff*make_rot(angl,centre);


angl=[0 0 params(3)].';
aff=aff*make_rot(angl,centre);


aff(1:3,4) = aff(1:3,4)+params(4:6).';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function n = norm2(v)
n=sqrt(sum(v.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function rot = make_rot(angl,centre)
rot=eye(4);
theta=norm2(angl);
if(theta < 1e-8)
    return;
end
axis=angl/theta;
x1=axis;
x2=[-axis(2) axis(1) 0].';
if(norm2(x2)<=0)
    x2 = [1 0 0].';
end
x2=x2/norm2(x2);
x3=cross(x1,x2);
x3=x3/norm2(x3);

b = [x2(:) x3(:) x1(:)];

rotcore=eye(3);
rotcore(1,1) = cos(theta);
rotcore(1,2) = sin(theta); 
rotcore(2,1) = -sin(theta);
rotcore(2,2) = cos(theta);

rot(1:3,1:3) = b*rotcore*(b.');
trans = (eye(3)-rot(1:3,1:3))*centre;
rot(1:3,4) = trans;

