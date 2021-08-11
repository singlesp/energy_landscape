function X = fillbbox(X,bbox,Xcropped)
sz=size(X);
bbox0=[ones(1,numel(sz)); sz];
bbox0(:,1:size(bbox,2))=bbox;

fnstr='';
for d = 1:size(bbox0,2)
    fnstr=[fnstr sprintf('bbox0(1,%d):bbox0(2,%d),',d,d)];
end
fnstr=sprintf('X(%s)=Xcropped;',fnstr(1:end-1));
eval(fnstr);
