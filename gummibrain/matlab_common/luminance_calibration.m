clear all;
close all;
clc;

numsteps = 6;
%numsteps = 20;
%numsteps = 40;

img0 = zeros(1,1,3);

figstat = figure('color',[0 0 0],'units','normalized');
txt = text(0,.5,'start','fontsize',20,'color',[1 1 1]);
axis off;

fig = figure('color',[0 0 0],'numbertitle','off','units','normalized');
imshow(img0);
set(gca,'position',[0 0 1 1]);
hold on;

levels = linspace(0,255,numsteps);

for i = 1:numel(levels)
    lev = levels(i);
    img = img0;
    %img(:,:,1) = lev/255; %red
    %img(:,:,2) = lev/255; %green
    %img(:,:,3) = lev/255; %blue
    img(:,:,:) = lev/255; %gray
    %img(:,:,[1 2]) = lev/255; %yellow
    %img(:,:,[1 3]) = lev/255; %magenta
    %img(:,:,[2 3]) = lev/255; %cyan
    
    %%%%%
    %img(:,:,1) = lev/255; %hue trace v=1
    %img(:,:,2) = 1;
    %img(:,:,3) = .7;
    %img = hsv2rgb(img);
    
    %%%%%
    %img(:,:,1) = lev/255; %hue trace v=.5
    %img(:,:,[2 3]) = .5;
    %img = hsv2rgb(img);
    
    set(0,'currentfigure',fig);
    
    imshow(img);
    s = sprintf('%d\n\n%.3f',i,lev/255);
    set(fig,'name',s);
    set(txt,'string',s);
    
    if(i < numel(levels))
        pause;
    end
end