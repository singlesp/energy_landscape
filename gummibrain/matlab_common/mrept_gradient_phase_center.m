% the positive directions are defined as
%      z
%     /
% x -
%    |
%    y
function [gx, gy, gz] = mrept_gradient_phase_center(A, d)

ndim=ndims(A);
size_input=size(A);
gx=zeros(size_input);
gy=zeros(size_input);
gz=zeros(size_input);

if ndim==2
    gx(:,2:size_input(2)-1)=A(:,1:size_input(2)-2)./A(:,3:size_input(2));
    gx(:,1)=gx(:,2);
    gx(:,size_input(2))=gx(:,size_input(2)-1);

    gy(2:size_input(1)-1,:)=A(3:size_input(1), :)./A(1:size_input(1)-2, :);
    gy(size_input(1),:)=gy(size_input(1)-1,:);
    gy(1,:)=gy(2,:);
end;

if ndim==3
    gx(:,2:size_input(2)-1,:)=A(:,1:size_input(2)-2, :)./A(:,3:size_input(2), :);
    gx(:,1,:)=gx(:,2,:);
    gx(:,size_input(2),:)=gx(:,size_input(2)-1,:);
    
    gy(2:size_input(1)-1,:,:)=A(3:size_input(1),:,:)./A(1:size_input(1)-2,:,:);
    gy(size_input(1),:,:)=gy(size_input(1)-1,:,:);
    gy(1,:,:)=gy(2,:,:);
        
    gz(:,:,2:size_input(3)-1)=A(:,:,3:size_input(3))./A(:,:,1:size_input(3)-2);
    gz(:,:,size_input(3))=gz(:,:,size_input(3)-1);
    gz(:,:,1)=gz(:,:,2);
end;

gx=angle(gx)/2/d(1);gy=angle(gy)/2/d(2);gz=angle(gz)/2/d(3);