function h = imshow_invisible(img)
sz=size(img);

%set(gcf,'visible','off','units','pixels','position',[0 0 max(sz(:))*2 max(sz())*2]);
set(gcf,'visible','off','units','pixels','position',[0 0 sz(2)+1 sz(1)+1]);
ax=axes('visible','off','units','pixels','position',[0 0 sz(2)+1 sz(1)+1]); 

hold on;


imgpad=zeros(size(img,1)+1,size(img,2)+1,size(img,3));
imgpad=cast(imgpad,'like',img);
imgpad(1:end-1,2:end,:)=img;
%imgpad(2:end-1,2:end-1,:)=img;


h=image(imgpad); 
axis tight image ij;
set(ax,'box','off','xtick',[],'ytick',[]);
