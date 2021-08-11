clear all;
close all;
clc;

%topo = imread('topomap.jpg');
topo = imread('flowerseg.jpg');
[w h tmp] = size(topo);

cr = reshape(topo(:,:,1),w*h,1);
cg = reshape(topo(:,:,2),w*h,1);
cb = reshape(topo(:,:,3),w*h,1);

%topocols = reshape(topo,w*h,3);
topocols = [cr cg cb];

plot3(cr,cg,cb,'.');

figure;

%imshow(topocols(k:k+200,:));
%axis off;

%subplot(3,1,1);
%plot(topocols(:,1));
%subplot(3,1,2);
%plot(topocols(:,2));
%subplot(3,1,3);
%plot(topocols(:,3));

k = (w*h- 10*h);
ind = k:k+h;

wh = 1:(w*h);

plot(wh(ind),topocols(ind,1),'r',wh(ind),topocols(ind,2),'g',wh(ind),topocols(ind,3),'b');

figure;

toposhow = zeros(w*h,1,3);
toposhow(:,1,:) = topocols;
imshow(toposhow(ind,:,:));

%%%%%%%%%%%%%%%%%%%%%%%%

%%
data = double(topocols);
L = distance(data');

return;
%%
d = data(topocols);

plot(d.X(:,1),d.X(:,2),'r.'); 
%d.Y=[];
pause;
[r,a]=train(spectral('sigma=5'),d);  %doesn't need classes
%[r a]=train(svm({kernel('rbf',1)}),d); %needs classes

I=find(r.X==1);clf;hold on;
plot(d.X(I,1),d.X(I,2),'r.');
I=find(r.X~=1);
plot(d.X(I,1),d.X(I,2),'b.');