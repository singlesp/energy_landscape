function h = stemoverlay(x,y,colors,varargin)

h = [];
[ymax,imax] = max(y,[],2);

if(isempty(colors))
    colors = get(gca,'colororder');
    colors = colors(1:size(y,2),:);
end

if(size(colors,1) == size(y,2)+1)
    cmix = colors(size(y,2)+1,:);
else
    cmix = mean(colors,1);
end

for i = 1:size(y,2)
    
    stem(x(imax==i),y(imax==i,i),'r',varargin{:},'color',colors(i,:));
    if(i==1)
        hold on;
    end  
end

for i = 1:size(y,2)
    stem(x(imax~=i),y(imax~=i,i),'m',varargin{:},'color',cmix);
end