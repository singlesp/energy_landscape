clear all;
close all;
clc;

sz = 100;
noise = .1;

p1 = rand(sz,3);
p1(:,3) = 0+noise*p1(:,3);
t1 = delaunay(p1(:,1),p1(:,2));

p2 = rand(sz,3);
p2(:,3) = 1+noise*p2(:,3);
t2 = delaunay(p2(:,1),p2(:,2));

trisurf(t1,p1(:,1),p1(:,2),p1(:,3));
hold on;
trisurf(t2,p2(:,1),p2(:,2),p2(:,3));

m = rand(sz,3);
m(:,3) = .5+noise*m(:,3);
plot3(m(:,1),m(:,2),m(:,3),'.');

%X = [p1; p2; m];
%n = size(p1,1)+size(p2,1);
X = [p1; m];
n = size(p1,1);

T = delaunayn(X);
%T = T - n;
%T = T(sum(T > 0,2) == 3,:);

%T = T(T > 0);
%T = reshape(T,numel(T)/3,3);

%[m,T]=meshcheckrepair(m,T);
%T=surfaceclean(T,m);

%figure;
%trisurf(T,m(:,1),m(:,2),m(:,3));

figure;
tetramesh(T,X);

%%
d = [-1 1];
[x,y,z] = meshgrid(d,d,d);  % A cube
x = [x(:);0];
y = [y(:);0];
z = [z(:);0];

% [x,y,z] are corners of a cube plus the center.
X = [x(:) y(:) z(:)];
X = X(1:end-1,:);
Tes = delaunayn(X,{'Qt','Qbb','Qc','Qz'});
tetramesh(Tes,X);
